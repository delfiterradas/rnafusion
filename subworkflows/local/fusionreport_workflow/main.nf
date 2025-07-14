include { FUSIONREPORT_DETECT      }     from '../../../modules/nf-core/fusionreport/detect/main'


workflow FUSIONREPORT_WORKFLOW {
    take:
        fusionreport_ref
        arriba_fusions
        starfusion_fusions
        fusioncatcher_fusions

    main:
        ch_versions = Channel.empty()
        ch_report = Channel.empty()
        ch_csv = Channel.empty()

        def ch_fusions = arriba_fusions
            .join(starfusion_fusions, failOnMismatch:true, failOnDuplicate:true)
            .join(fusioncatcher_fusions, failOnMismatch:true, failOnDuplicate:true)

        FUSIONREPORT_DETECT(ch_fusions, fusionreport_ref, params.tools_cutoff)
        ch_fusion_list = FUSIONREPORT_DETECT.out.fusion_list
        ch_fusion_list_filtered = FUSIONREPORT_DETECT.out.fusion_list_filtered
        ch_versions = ch_versions.mix(FUSIONREPORT_DETECT.out.versions)
        ch_report = FUSIONREPORT_DETECT.out.report
        ch_csv = FUSIONREPORT_DETECT.out.csv

    emit:
        versions                 = ch_versions
        fusion_list              = ch_fusion_list
        fusion_list_filtered     = ch_fusion_list_filtered
        report                   = ch_report
        csv                      = ch_csv

}
