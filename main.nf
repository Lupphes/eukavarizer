#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/eukavarizer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/eukavarizer
    Website: https://nf-co.re/eukavarizer
    Slack  : https://nfcore.slack.com/channels/eukavarizer
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SEQRETRIEVAL            } from './workflows/seqretrieval'
include { EUKAVARIZER             } from './workflows/eukavarizer'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_eukavarizer_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_eukavarizer_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Define the workflow
workflow NFCORE_EUKAVARIZER {
    take:
        ch_taxonomy_id
        ch_outdir
        ch_local_sequences_dir
        ch_local_refseq_path

    main:
        SEQRETRIEVAL (
            ch_taxonomy_id,
            ch_outdir,
            ch_local_sequences_dir,
            ch_local_refseq_path
        )

        //
        // WORKFLOW: Run pipeline
        //
        // EUKAVARIZER (
        //     SEQRETRIEVAL.out.fastq_files,
        //     SEQRETRIEVAL.out.refseq_path,
        //     ch_local_sequences_dir,
        //     ch_local_refseq_path
        // )

    emit:
        multiqc_report = "test" //EUKAVARIZER.out.multiqc_report
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    //
    // Input channels for the workflow
    //
    main:
    ch_taxonomy_id         = params.taxonomy_id ? Channel.value(params.taxonomy_id) : Channel.empty()
    outdir_channel         = params.outdir ? Channel.value(params.outdir) : Channel.empty()
    ch_local_sequences_dir = params.local_sequences_dir ? Channel.value(params.local_sequences_dir) : Channel.empty()
    ch_local_refseq_path   = params.local_refseq_path ? Channel.value(params.local_refseq_path) : Channel.empty()


    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        outdir_channel
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_EUKAVARIZER (
        ch_taxonomy_id,
        outdir_channel,
        ch_local_sequences_dir,
        ch_local_refseq_path
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
        NFCORE_EUKAVARIZER.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
