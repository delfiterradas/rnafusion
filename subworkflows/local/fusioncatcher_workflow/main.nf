include { FUSIONCATCHER } from '../../../modules/local/fusioncatcher/detect/main'

workflow FUSIONCATCHER_WORKFLOW {
    take:
        reads                   // channel [ meta, [ fastqs ] ]
        fusioncatcher_ref       // channel [ meta, path       ]
        run_fusioncatcher       // boolean
        all                     // boolean
        fusioninspector_only    // boolean
        fusioncatcher_fusions   // path, string

    main:
        ch_versions   = Channel.empty()
        ch_dummy_file = file("$baseDir/assets/dummy_file_fusioncatcher.txt", checkIfExists: true)

        if (( run_fusioncatcher || all) && !fusioninspector_only ) {
            if (fusioncatcher_fusions){

                ch_fusioncatcher_fusions = reads.combine(Channel.value(file(fusioncatcher_fusions, checkIfExists:true)))
                                            .map { meta, reads, fusions -> [ meta, fusions ] }
            } else {

                FUSIONCATCHER (
                    reads,
                    fusioncatcher_ref
                )
                ch_fusioncatcher_fusions = FUSIONCATCHER.out.fusions
                ch_versions = ch_versions.mix(FUSIONCATCHER.out.versions)
            }
        }
        else {
            ch_fusioncatcher_fusions = reads.combine(Channel.value(file(ch_dummy_file, checkIfExists:true)))
                                        .map { meta, reads, fusions -> [ meta, fusions ] }
        }

    emit:
        fusions  = ch_fusioncatcher_fusions
        versions = ch_versions
    }

