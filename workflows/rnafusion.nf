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

    main:

    def ch_versions = Channel.empty()
    def ch_multiqc_files = Channel.empty()

    //
    // Create references if necessary
    //

    BUILD_REFERENCES()
    ch_versions = ch_versions.mix(BUILD_REFERENCES.out.versions)

    if (!params.references_only) { // TODO: Remove this temporary parameter when we have a full-working GitHub nf-test

        //
        // QC from FASTQ files
        //
        FASTQC (
            ch_samplesheet
        )
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
        ch_versions = ch_versions.mix(FASTQC.out.versions)


        //
        // Trimming
        //
        TRIM_WORKFLOW (
            ch_samplesheet,
            Channel.value(params.adapter_fasta),
            params.fastp_trim
        )
        ch_reads = TRIM_WORKFLOW.out.ch_reads_all
        ch_versions = ch_versions.mix(TRIM_WORKFLOW.out.versions)

        SALMON_QUANT( ch_reads, BUILD_REFERENCES.out.ch_salmon_index.map{ it -> it[1] }, BUILD_REFERENCES.out.ch_gtf.map{ it -> it[1] }, [], false, 'A')
        ch_multiqc_files = ch_multiqc_files.mix(SALMON_QUANT.out.json_info.collect{it[1]})
        ch_versions = ch_versions.mix(SALMON_QUANT.out.versions)

        //
        // Run STAR alignment
        //

        def ch_aligned_reads = Channel.empty()
        def ch_star_junctions = Channel.empty()
        def ch_star_split_junctions = Channel.empty()
        if((params.arriba || params.all || params.ctatsplicing || params.starfusion) && !params.fusioninspector_only) {
            FASTQ_ALIGN_STAR(
                ch_reads,
                BUILD_REFERENCES.out.ch_starindex_ref,
                BUILD_REFERENCES.out.ch_gtf,
                BUILD_REFERENCES.out.ch_fasta,
                BUILD_REFERENCES.out.ch_fai,
                params.star_ignore_sjdbgtf,
                params.seq_platform,
                params.seq_center,
                params.cram
            )
            ch_versions = ch_versions.mix(FASTQ_ALIGN_STAR.out.versions)
            ch_aligned_reads = FASTQ_ALIGN_STAR.out.bam_bai
            ch_star_junctions = FASTQ_ALIGN_STAR.out.junctions
            ch_star_split_junctions = FASTQ_ALIGN_STAR.out.spl_junc_tabs
        }

        //
        // Run CTAT-SPLICING
        //

        if((params.ctatsplicing || params.all)) {
            CTATSPLICING_WORKFLOW(
                ch_star_split_junctions,
                ch_star_junctions,
                ch_aligned_reads.map { meta, bam, _bai -> [ meta, bam ]},
                BUILD_REFERENCES.out.ch_starfusion_ref
            )
            ch_versions = ch_versions.mix(CTATSPLICING_WORKFLOW.out.versions)
        }

        //
        // SUBWORKFLOW:  Run STAR alignment and Arriba
        //

        // TODO: add params.seq_platform and pass it as argument to arriba_workflow
        // TODO: improve how params.arriba_fusions would avoid running arriba module. Maybe imputed from samplesheet?
        // TODO: same as above, but with ch_arriba_fusion_fail. It's currently replaces by a dummy file

        ARRIBA_WORKFLOW (
            ch_aligned_reads.map { meta, bam, _bai -> [ meta, bam ]},
            BUILD_REFERENCES.out.ch_gtf,
            BUILD_REFERENCES.out.ch_fasta,
            BUILD_REFERENCES.out.ch_arriba_ref_blacklist,
            BUILD_REFERENCES.out.ch_arriba_ref_cytobands,
            BUILD_REFERENCES.out.ch_arriba_ref_known_fusions,
            BUILD_REFERENCES.out.ch_arriba_ref_protein_domains,
            params.arriba,                   // boolean
            params.all,                      // boolean
            params.fusioninspector_only,     // boolean
            params.arriba_fusions,           // path
        )
        ch_versions = ch_versions.mix(ARRIBA_WORKFLOW.out.versions)


        //Run STAR fusion
        STARFUSION_WORKFLOW (
            ch_aligned_reads.map { meta, bam, _bai -> [ meta, bam ]},
            ch_star_junctions,
            BUILD_REFERENCES.out.ch_starfusion_ref
        )
        ch_versions = ch_versions.mix(STARFUSION_WORKFLOW.out.versions)

        //Run fusioncatcher
        FUSIONCATCHER_WORKFLOW (
            ch_reads,
            BUILD_REFERENCES.out.ch_fusioncatcher_ref,       // channel [ meta, path       ]
            params.run_fusioncatcher,
            params.all,
            params.fusioninspector_only,
            params.fusioncatcher_fusions
        )
        ch_versions = ch_versions.mix(FUSIONCATCHER_WORKFLOW.out.versions)


        //Run stringtie
        STRINGTIE_WORKFLOW (
            ch_aligned_reads.map { meta, bam, _bai -> [meta, bam]},
            BUILD_REFERENCES.out.ch_gtf
        )
        ch_versions = ch_versions.mix(STRINGTIE_WORKFLOW.out.versions)


        //Run fusion-report
        FUSIONREPORT_WORKFLOW (
            ch_reads,
            BUILD_REFERENCES.out.ch_fusionreport_ref,
            ARRIBA_WORKFLOW.out.fusions,
            STARFUSION_WORKFLOW.out.fusions,
            FUSIONCATCHER_WORKFLOW.out.fusions
        )
        ch_versions = ch_versions.mix(FUSIONREPORT_WORKFLOW.out.versions)

        //Run fusionInpector
        FUSIONINSPECTOR_WORKFLOW (
            ch_reads,
            FUSIONREPORT_WORKFLOW.out.fusion_list,
            FUSIONREPORT_WORKFLOW.out.fusion_list_filtered,
            FUSIONREPORT_WORKFLOW.out.report,
            FUSIONREPORT_WORKFLOW.out.csv,
            STARFUSION_WORKFLOW.out.ch_bam_sorted_indexed,
            BUILD_REFERENCES.out.ch_gtf,
            BUILD_REFERENCES.out.ch_arriba_ref_protein_domains,
            BUILD_REFERENCES.out.ch_arriba_ref_cytobands,
            BUILD_REFERENCES.out.ch_hgnc_ref,
            BUILD_REFERENCES.out.ch_hgnc_date
        )
        ch_versions = ch_versions.mix(FUSIONINSPECTOR_WORKFLOW.out.versions)


        //QC
        QC_WORKFLOW (
            ch_aligned_reads,
            BUILD_REFERENCES.out.ch_refflat,
            BUILD_REFERENCES.out.ch_fasta,
            BUILD_REFERENCES.out.ch_fai,
            BUILD_REFERENCES.out.ch_rrna_interval
        )
        ch_versions = ch_versions.mix(QC_WORKFLOW.out.versions)

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

    if (!params.references_only) { // TODO: Remove this temporary parameter when we have a full-working GitHub nf-test
        ch_multiqc_files                      = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files                      = ch_multiqc_files.mix(TRIM_WORKFLOW.out.ch_fastp_html.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files                      = ch_multiqc_files.mix(TRIM_WORKFLOW.out.ch_fastp_json.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files                      = ch_multiqc_files.mix(TRIM_WORKFLOW.out.ch_fastqc_trimmed.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files                      = ch_multiqc_files.mix(FASTQ_ALIGN_STAR.out.log_final.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files                      = ch_multiqc_files.mix(FASTQ_ALIGN_STAR.out.gene_count.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files                      = ch_multiqc_files.mix(QC_WORKFLOW.out.rnaseq_metrics.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files                      = ch_multiqc_files.mix(QC_WORKFLOW.out.duplicate_metrics.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files                      = ch_multiqc_files.mix(QC_WORKFLOW.out.insertsize_metrics.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files                      = ch_multiqc_files.mix(FUSIONINSPECTOR_WORKFLOW.out.ch_arriba_visualisation.collect{it[1]}.ifEmpty([]))
    }

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )



    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
