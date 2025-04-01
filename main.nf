#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/rnafusion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/rnafusion
    Website: https://nf-co.re/rnafusion
    Slack  : https://nfcore.slack.com/channels/rnafusion
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_rnafusion_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_rnafusion_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_rnafusion_pipeline'
include { RNAFUSION               } from './workflows/rnafusion'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        args,
        params.outdir,
    )

    // A list of all tools available in the pipeline. Add tools here when a new tool has been added
    def all_tools = [
        "arriba",
        "ctatsplicing",
        "fusioncatcher",
        "starindex",
        "starfusion",
        "stringtie",
        "fusionreport",
        "fastp",
        "salmon",
        "fusioninspector"
    ]

    def tools = params.all ? all_tools : all_tools.findAll { tool ->
        params.get(tool)
    }

    //
    // WORKFLOW: Run main workflow
    //
    RNAFUSION(
        PIPELINE_INITIALISATION.out.samplesheet,
        tools
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        RNAFUSION.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
