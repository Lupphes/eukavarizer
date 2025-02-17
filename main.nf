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
        ch_sequence_dir
        ch_genome_file

    main:

        println "Running nf-core/eukavarizer pipeline"

        //
        // STEP 1: Retrieve Sequences & Reference Genome
        //
        seqretrieval_results = SEQRETRIEVAL(
            ch_taxonomy_id,
            ch_outdir,
            ch_sequence_dir,
            ch_genome_file
        )

        //
        // STEP 2: Run Main Analysis Pipeline
        //
        eukavarizer_results = EUKAVARIZER(
            seqretrieval_results.genome_file,
            seqretrieval_results.fastq_files
        )

    emit:
        multiqc_report = eukavarizer_results.multiqc_report
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
        ch_taxonomy_id         = Channel.value(params.taxonomy_id ?: '4932')
        ch_outdir              = Channel.value(params.outdir ?: './output')
        sequence_dir           = Channel.value(params.sequence_dir ?: "")
        ch_genome_file         = Channel.value(params.genome_file ?: "")

        //
        // SUBWORKFLOW: Run initialisation tasks
        //
        PIPELINE_INITIALISATION (
            params.version,
            params.validate_params,
            params.monochrome_logs,
            args,
            ch_outdir
        )

        //
        // WORKFLOW: Run main workflow
        //
        NFCORE_EUKAVARIZER (
            ch_taxonomy_id,
            ch_outdir,
            sequence_dir,
            ch_genome_file
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
