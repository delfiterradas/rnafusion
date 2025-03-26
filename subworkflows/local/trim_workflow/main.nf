

include { FASTP                      } from '../../../modules/nf-core/fastp/main'
include { FASTQC as FASTQC_FOR_FASTP } from '../../../modules/nf-core/fastqc/main'

workflow TRIM_WORKFLOW {

    take:
        reads           // channel [ meta, [ fastq files ] ]
        adapter_fasta   // channel [ path ]
        fastp_trim      // boolean

    main:
        ch_versions       = Channel.empty()
        ch_fastp_html     = Channel.empty()
        ch_fastp_json     = Channel.empty()
        ch_fastqc_trimmed = Channel.empty()

        if ( fastp_trim ) {
            FASTP(reads, adapter_fasta.ifEmpty( [] ), false, false, false)
            ch_versions = ch_versions.mix(FASTP.out.versions)

            FASTQC_FOR_FASTP(FASTP.out.reads)
            ch_versions = ch_versions.mix(FASTQC_FOR_FASTP.out.versions)

            ch_reads_all           = FASTP.out.reads
            ch_reads_fusioncatcher = ch_reads_all
            ch_fastp_html          = FASTP.out.html
            ch_fastp_json          = FASTP.out.json
            ch_fastqc_trimmed      = FASTQC_FOR_FASTP.out.zip

        }
        else {
            ch_reads_all           = reads
            ch_reads_fusioncatcher = reads
        }

    emit:
        ch_reads_all            // Channel [ meta, [reads]   ]
        ch_reads_fusioncatcher  // Channel [ meta, [reads]   ]
        ch_fastp_html           // Channel [ meta, path_html ]
        ch_fastp_json           // Channel [ meta, path_json ]
        ch_fastqc_trimmed       // Channel [ meta, path_zip  ]
        versions = ch_versions  // Channel [ versions        ]
    }

