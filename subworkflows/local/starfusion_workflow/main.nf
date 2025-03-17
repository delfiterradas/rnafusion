include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_FOR_STARFUSION       }   from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_FOR_STARFUSION_CRAM  }   from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_FOR_STARFUSION         }   from '../../../modules/nf-core/samtools/view/main'
include { STARFUSION                                            }   from '../../../modules/local/starfusion/detect/main'
include { CTATSPLICING_WORKFLOW                                 }   from '../ctatsplicing_workflow'

workflow STARFUSION_WORKFLOW {
    take:
        reads               // channel: [ meta, bam ]
        junctions           // channel: [ meta, junctions ]
        ch_starfusion_ref

    main:
        def ch_versions = Channel.empty()
        def ch_align = Channel.empty()
        def ch_starfusion_fusions = Channel.empty()
        def bam_sorted_indexed = Channel.empty()

        ch_dummy_file = file("$baseDir/assets/dummy_file_starfusion.txt", checkIfExists: true)

        if ((params.starfusion || params.all || params.stringtie) && !params.fusioninspector_only) {
            if (params.starfusion_fusions){
                ch_starfusion_fusions = reads.combine(Channel.value(file(params.starfusion_fusions, checkIfExists:true)))
                                        .map { it -> [ it[0], it[2] ] }
            } else {
                reads_junction = reads.join(junctions)  // TODO: This join is not needed as STARFUSION can simply read from the junction file: https://github.com/STAR-Fusion/STAR-Fusion/wiki#alternatively-kickstart-mode-running-star-yourself-and-then-running-star-fusion-using-the-existing-outputs

                if (params.starfusion || params.all){
                    STARFUSION( reads_junction, ch_starfusion_ref.map { it -> it[1] })
                    ch_versions = ch_versions.mix(STARFUSION.out.versions)
                    ch_starfusion_fusions = STARFUSION.out.fusions
                }
            }
        }
        else {
            ch_starfusion_fusions = reads.combine(Channel.value(file(ch_dummy_file, checkIfExists:true)))
                                    .map { it -> [ it[0], it[2] ] }
        }
    emit:
        fusions               = ch_starfusion_fusions
        ch_bam_sorted         = ch_align
        ch_bam_sorted_indexed = bam_sorted_indexed
        versions              = ch_versions
    }

