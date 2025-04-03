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
        starfusion_fusions

    main:
        def ch_versions = Channel.empty()
        def ch_starfusion_fusions = Channel.empty()

        if (starfusion_fusions){
            ch_starfusion_fusions = reads.combine(Channel.value(file(starfusion_fusions, checkIfExists:true)))
                                    .map { it -> [ it[0], it[2] ] }
        } else {
            reads_junction = reads.join(junctions)  // TODO: This join is not needed as STARFUSION can simply read from the junction file: https://github.com/STAR-Fusion/STAR-Fusion/wiki#alternatively-kickstart-mode-running-star-yourself-and-then-running-star-fusion-using-the-existing-outputs

            STARFUSION( reads_junction, ch_starfusion_ref.map { it -> it[1] })
            ch_versions = ch_versions.mix(STARFUSION.out.versions)
            ch_starfusion_fusions = STARFUSION.out.fusions

        }
    emit:
        fusions               = ch_starfusion_fusions
        versions              = ch_versions
    }

