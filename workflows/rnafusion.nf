/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BUILD_REFERENCES              }   from '../subworkflows/local/build_references'
include { CAT_FASTQ                     }   from '../modules/nf-core/cat/fastq/main'
include { TRIM_WORKFLOW                 }   from '../subworkflows/local/trim_workflow/main'
include { ARRIBA_WORKFLOW               }   from '../subworkflows/local/arriba_workflow'
include { QC_WORKFLOW                   }   from '../subworkflows/local/qc_workflow'
include { STARFUSION_WORKFLOW           }   from '../subworkflows/local/starfusion_workflow'
include { STRINGTIE_WORKFLOW            }   from '../subworkflows/local/stringtie_workflow/main'
include { FUSIONCATCHER_WORKFLOW        }   from '../subworkflows/local/fusioncatcher_workflow'
include { FUSIONINSPECTOR_WORKFLOW      }   from '../subworkflows/local/fusioninspector_workflow'
include { FUSIONREPORT_WORKFLOW         }   from '../subworkflows/local/fusionreport_workflow'
include { FASTQC                        }   from '../modules/nf-core/fastqc/main'
include { MULTIQC                       }   from '../modules/nf-core/multiqc/main'
include { STAR_ALIGN                    }   from '../modules/nf-core/star/align/main'
include { SALMON_QUANT                  }   from '../modules/nf-core/salmon/quant/main'
include { paramsSummaryMap              }   from 'plugin/nf-schema'
include { FASTQ_ALIGN_STAR              }   from '../subworkflows/local/fastq_align_star'
include { CTATSPLICING_WORKFLOW         }   from '../subworkflows/local/ctatsplicing_workflow'
include { paramsSummaryMultiqc          }   from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML        }   from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText        }   from '../subworkflows/local/utils_nfcore_rnafusion_pipeline'
include { validateInputSamplesheet      }   from '../subworkflows/local/utils_nfcore_rnafusion_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RNAFUSION {


    take:
    ch_samplesheet // channel: samplesheet read in from --input
    tools          // list: a list of tools to run

    main:

    def ch_versions = Channel.empty()
    def ch_multiqc_files = Channel.empty()

    //
    // Create references if necessary
    //

    BUILD_REFERENCES(tools)
    ch_versions = ch_versions.mix(BUILD_REFERENCES.out.versions)

    if (!params.references_only) {

        //
        // QC from FASTQ files
        //

        if(!params.skip_qc) {
            FASTQC (
                ch_samplesheet
            )
            ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
            ch_versions = ch_versions.mix(FASTQC.out.versions)
        }

        //
        // SUBWORKFLOW: Trimming
        //

        def ch_reads = Channel.empty()
        if(tools.contains("fastp")) {
            def ch_adapter_fasta = params.adapter_fasta ? Channel.fromPath(params.adapter_fasta).collect() : []
            TRIM_WORKFLOW (
                ch_samplesheet,
                ch_adapter_fasta,
                params.skip_qc
            )
            ch_reads = TRIM_WORKFLOW.out.ch_reads_all
            ch_versions      = ch_versions.mix(TRIM_WORKFLOW.out.versions)
            ch_multiqc_files = ch_multiqc_files.mix(TRIM_WORKFLOW.out.ch_fastp_html.collect{it[1]})
            ch_multiqc_files = ch_multiqc_files.mix(TRIM_WORKFLOW.out.ch_fastp_json.collect{it[1]})
            ch_multiqc_files = ch_multiqc_files.mix(TRIM_WORKFLOW.out.ch_fastqc_trimmed.collect{it[1]})
        } else {
            ch_reads = ch_samplesheet
        }

        if(tools.contains("salmon") && !params.skip_qc) {
            SALMON_QUANT(
                ch_reads,
                BUILD_REFERENCES.out.salmon_index,
                BUILD_REFERENCES.out.gtf.map{ it -> it[1] },
                [],
                false,
                'A'
            )
            ch_multiqc_files = ch_multiqc_files.mix(SALMON_QUANT.out.json_info.collect{it[1]})
            ch_versions      = ch_versions.mix(SALMON_QUANT.out.versions)
        }
        //
        // Run STAR alignment
        //

        def ch_aligned_reads = Channel.empty()
        def ch_star_junctions = Channel.empty()
        def ch_star_split_junctions = Channel.empty()
        if(tools.intersect(["ctatsplicing", "arriba", "starfusion", "stringtie"])) {
            FASTQ_ALIGN_STAR(
                ch_reads,
                BUILD_REFERENCES.out.starindex_ref,
                BUILD_REFERENCES.out.gtf,
                BUILD_REFERENCES.out.fasta,
                BUILD_REFERENCES.out.fai,
                params.star_ignore_sjdbgtf,
                params.seq_platform,
                params.seq_center,
                params.cram
            )
            ch_versions = ch_versions.mix(FASTQ_ALIGN_STAR.out.versions)
            ch_aligned_reads = FASTQ_ALIGN_STAR.out.bam_bai
            ch_star_junctions = FASTQ_ALIGN_STAR.out.junctions
            ch_star_split_junctions = FASTQ_ALIGN_STAR.out.spl_junc_tabs
            ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_STAR.out.log_final.collect{it[1]}.ifEmpty([]))
            ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_STAR.out.gene_count.collect{it[1]}.ifEmpty([]))
        }

        //
        // Run CTAT-SPLICING
        //

        if(tools.contains("ctatsplicing")) {
            CTATSPLICING_WORKFLOW(
                ch_star_split_junctions,
                ch_star_junctions,
                ch_aligned_reads,
                BUILD_REFERENCES.out.starfusion_ref
            )
            ch_versions = ch_versions.mix(CTATSPLICING_WORKFLOW.out.versions)
        }

        //
        // SUBWORKFLOW: Run STAR alignment and Arriba
        //

        // TODO: improve how params.arriba_fusions would avoid running arriba module. Maybe imputed from samplesheet?
        // TODO: same as above, but with ch_arriba_fusion_fail. It's currently replaces by a dummy file

        def fusions_created = false
        def ch_arriba_fusions = ch_reads.map { it -> [it[0], []] } // Set arriba fusions to empty by default
        if(tools.contains("arriba")) {
            fusions_created = true
            ARRIBA_WORKFLOW (
                ch_aligned_reads.map { meta, bam, _bai -> [ meta, bam ]},
                BUILD_REFERENCES.out.gtf,
                BUILD_REFERENCES.out.fasta,
                BUILD_REFERENCES.out.arriba_ref_blacklist,
                BUILD_REFERENCES.out.arriba_ref_cytobands,
                BUILD_REFERENCES.out.arriba_ref_known_fusions,
                BUILD_REFERENCES.out.arriba_ref_protein_domains,
                params.arriba_fusions
            )
            ch_versions = ch_versions.mix(ARRIBA_WORKFLOW.out.versions)
            ch_arriba_fusions = ARRIBA_WORKFLOW.out.fusions
        }

        //
        // SUBWORKFLOW: Run STAR alignment and StarFusion
        //

        def ch_starfusion_fusions = ch_reads.map { it -> [it[0], []] } // Set starfusion fusions to empty by default
        if(tools.contains("starfusion")) {
            fusions_created = true
            STARFUSION_WORKFLOW (
                ch_aligned_reads.map { meta, bam, _bai -> [ meta, bam ]},
                ch_star_junctions,
                BUILD_REFERENCES.out.starfusion_ref,
                params.starfusion_fusions
            )
            ch_versions             = ch_versions.mix(STARFUSION_WORKFLOW.out.versions)
            ch_starfusion_fusions   = STARFUSION_WORKFLOW.out.fusions
        }

        //
        // SUBWORKFLOW: Run FusionCatcher
        //

        def ch_fusioncatcher_fusions = ch_reads.map { it -> [it[0], []] } // Set fusioncatcher fusions to empty by default
        if(tools.contains("fusioncatcher")) {
            fusions_created = true
            FUSIONCATCHER_WORKFLOW (
                ch_reads,
                BUILD_REFERENCES.out.fusioncatcher_ref,       // channel [ meta, path       ]
                params.fusioncatcher_fusions
            )
            ch_versions = ch_versions.mix(FUSIONCATCHER_WORKFLOW.out.versions)
        }

        //
        // SUBWORKFLOW: Run Stringtie
        //

        if(tools.contains("stringtie")) {
            STRINGTIE_WORKFLOW (
                ch_aligned_reads.map { meta, bam, _bai -> [meta, bam]},
                BUILD_REFERENCES.out.gtf
            )
            ch_versions = ch_versions.mix(STRINGTIE_WORKFLOW.out.versions)
        }

        //
        // SUBWORKFLOW: Run FusionReport
        //

        def ch_fusion_list = Channel.empty()
        def ch_fusion_list_filtered = Channel.empty()
        def ch_fusionreport_report = Channel.empty()
        def ch_fusionreport_csv = Channel.empty()
        if(!params.skip_vis && tools.contains("fusionreport")) {
            if (!fusions_created) {
                error("Could not find any fusion files. Please generate some with --arriba, --starfusion and/or --fusioncatcher")
            }
            FUSIONREPORT_WORKFLOW (
                ch_reads,
                BUILD_REFERENCES.out.fusionreport_ref,
                ch_arriba_fusions,
                ch_starfusion_fusions,
                ch_fusioncatcher_fusions
            )
            ch_versions             = ch_versions.mix(FUSIONREPORT_WORKFLOW.out.versions)
            ch_fusion_list          = FUSIONREPORT_WORKFLOW.out.fusion_list
            ch_fusion_list_filtered = FUSIONREPORT_WORKFLOW.out.fusion_list_filtered
            ch_fusionreport_report  = FUSIONREPORT_WORKFLOW.out.report
            ch_fusionreport_csv     = FUSIONREPORT_WORKFLOW.out.csv
        } else if(params.fusioninspector_fusions) {
            def input_fusions       = file(params.fusioninspector_fusions, checkIfExists:true)
            ch_fusion_list          = ch_reads.map { it -> [ it[0], input_fusions ] }
            ch_fusion_list_filtered = ch_fusion_list
            ch_fusionreport_csv     = null
            ch_fusionreport_report  = null
        } else if(tools.contains("fusioninspector")) {
            error("Could not find any valid fusions for fusioninspector input. Please provide some via --fusioninspector_fusions or generate them with --arriba, --starfusion and/or --fusioncatcher with --skip_vis disabled and --fusionreport enabled")
        }

        //
        // SUBWORKFLOW: Run FusionInspector
        //

        if (tools.contains("fusioninspector")) {
            FUSIONINSPECTOR_WORKFLOW (
                ch_reads,
                ch_fusion_list,
                ch_fusion_list_filtered,
                ch_fusionreport_report,
                ch_fusionreport_csv,
                ch_aligned_reads,
                BUILD_REFERENCES.out.gtf,
                BUILD_REFERENCES.out.arriba_ref_protein_domains,
                BUILD_REFERENCES.out.arriba_ref_cytobands,
                BUILD_REFERENCES.out.hgnc_ref,
                BUILD_REFERENCES.out.hgnc_date,
                BUILD_REFERENCES.out.starfusion_ref,
                params.skip_vis,
                params.skip_vcf
            )
            ch_versions      = ch_versions.mix(FUSIONINSPECTOR_WORKFLOW.out.versions)
            ch_multiqc_files = ch_multiqc_files.mix(FUSIONINSPECTOR_WORKFLOW.out.ch_arriba_visualisation.collect{it[1]}.ifEmpty([]))
        }

        //
        // SUBWORKFLOW: Run QC
        //

        if(!params.skip_qc) {
            QC_WORKFLOW (
                ch_aligned_reads.map { meta, bam, _bai -> [meta, bam] },
                BUILD_REFERENCES.out.refflat,
                BUILD_REFERENCES.out.fasta,
                BUILD_REFERENCES.out.fai,
                BUILD_REFERENCES.out.rrna_interval
            )
            ch_versions      = ch_versions.mix(QC_WORKFLOW.out.versions)
            ch_multiqc_files = ch_multiqc_files.mix(QC_WORKFLOW.out.rnaseq_metrics.collect{it[1]})
            ch_multiqc_files = ch_multiqc_files.mix(QC_WORKFLOW.out.duplicate_metrics.collect{it[1]})
            ch_multiqc_files = ch_multiqc_files.mix(QC_WORKFLOW.out.insertsize_metrics.collect{it[1]})
        }
    }
    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'rnafusion_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //

    def ch_multiqc_output = Channel.empty()
    if(!params.skip_qc && !params.references_only) {
        ch_multiqc_config        = Channel.fromPath(
            "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
        ch_multiqc_custom_config = params.multiqc_config ?
            Channel.fromPath(params.multiqc_config, checkIfExists: true) :
            Channel.empty()
        ch_multiqc_logo          = params.multiqc_logo ?
            Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
            Channel.empty()

        summary_params      = paramsSummaryMap(
            workflow, parameters_schema: "nextflow_schema.json")
        ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
        ch_multiqc_files = ch_multiqc_files.mix(
            ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
            file(params.multiqc_methods_description, checkIfExists: true) :
            file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
        ch_methods_description                = Channel.value(
            methodsDescriptionText(ch_multiqc_custom_methods_description))

        ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
        ch_multiqc_files = ch_multiqc_files.mix(
            ch_methods_description.collectFile(
                name: 'methods_description_mqc.yaml',
                sort: true
            )
        )

        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList(),
            [],
            []
        )
        ch_multiqc_output = MULTIQC.out.report.toList()
    }


    emit:
    multiqc_report =  ch_multiqc_output // channel: /path/to/multiqc_report.html
    versions       = ch_versions        // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
