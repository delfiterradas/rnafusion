include { ARRIBA_ARRIBA                               } from '../../../modules/nf-core/arriba/arriba/main'

workflow ARRIBA_WORKFLOW {
    take:
        reads                           // channel [ meta, [ bam ]            ]
        ch_gtf                          // channel [ meta, path_gtf           ]
        ch_fasta                        // channel [ meta, path_fasta         ]
        ch_arriba_ref_blacklist         // channel [ meta, path_blacklist     ]
        ch_arriba_ref_cytobands         // channel [ meta, path_cytobands     ]
        ch_arriba_ref_known_fusions     // channel [ meta, path_known_fusions ]
        ch_arriba_ref_protein_domains   // channel [ meta, path_proteins      ]
        arriba                          // boolean
        all                             // boolean
        fusioninspector_only            // boolean
        arriba_fusions                  // path

    main:

        def ch_versions   = Channel.empty()
        def ch_dummy_file = file("$projectDir/assets/dummy_file_arriba.txt", checkIfExists: true)

        if (( arriba || all ) && !fusioninspector_only) {

            if ( arriba_fusions ) {

                ch_arriba_fusions = reads.combine( Channel.value( file( arriba_fusions, checkIfExists: true ) ) )
                    .map { it -> [ it[0], it[2] ] }
                ch_arriba_fusion_fail = ch_dummy_file

            } else {

                ARRIBA_ARRIBA (
                    reads,
                    ch_fasta,
                    ch_gtf,
                    ch_arriba_ref_blacklist,
                    ch_arriba_ref_known_fusions,
                    ch_arriba_ref_cytobands,
                    ch_arriba_ref_protein_domains
                )

                ch_versions = ch_versions.mix(ARRIBA_ARRIBA.out.versions)

                ch_arriba_fusions     = ARRIBA_ARRIBA.out.fusions
                ch_arriba_fusion_fail = ARRIBA_ARRIBA.out.fusions_fail.map{ it -> return it[1] }
            }

        } else {

            ch_arriba_fusions = reads
                .combine(Channel.value( file(ch_dummy_file, checkIfExists: true ) ) )
                .map { it -> [ it[0], it[2] ] }

            ch_arriba_fusion_fail = ch_dummy_file
        }

    emit:
        fusions      = ch_arriba_fusions        // channel [ meta, path_fusions   ]
        fusions_fail = ch_arriba_fusion_fail    // channel [ path, fusions_failed ]
        versions     = ch_versions              // channel [ versions             ]
    }

