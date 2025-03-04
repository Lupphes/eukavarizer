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
include { ANALYSIS_FAST_MULTI     } from './workflows/analysis_fast_multi'
include { INPUT_GENERATION        } from './workflows/input_generation'
include { REPORT_GENERATION       } from './workflows/report_generation'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_eukavarizer_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_eukavarizer_pipeline'
include { STRUCTURAL_VARIANT_CALLING } from './workflows/structural_variant_calling'


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
        ch_sequence_dir_abs
        ch_sequence_dir
        ch_genome_file

    main:

        view("🚀 Running nf-core/eukavarizer pipeline")
        //
        // STEP 1: Retrieve Sequences & Reference Genome
        //
        seqretrieval_results = SEQRETRIEVAL(
            ch_taxonomy_id,
            ch_outdir,
            ch_sequence_dir_abs,
            ch_genome_file
        )

        seqretrieval_results.grouped_fastqs.view { "DEBUG: seqretrieval_results.grouped_fastqs -> ${it}" }

        //
        // STEP 2: Run Main Analysis Pipeline
        //
        data_analysis_results = ANALYSIS_FAST_MULTI(
            seqretrieval_results.genome_file,
            seqretrieval_results.fastq_files
        )

        // //
        // // STEP 3: Preprocess and generate INPUT for Structural Analysis
        // // with alignments
        // //
        view("🔹 Checking inputs for INPUT_GENERATION")

        seqretrieval_results.fastq_files.view { "DEBUG: seqretrieval_results.fastq_files -> ${it}" }
        seqretrieval_results.genome_file.view { "DEBUG: seqretrieval_results.genome_file -> ${it}" }

        input_generation_results = INPUT_GENERATION(
            seqretrieval_results.grouped_fastqs,
            seqretrieval_results.genome_file
        )


        // // Debugging: Output BAM files before variant calling
        input_generation_results.bam_files.view { "DEBUG BAM FILES BEFORE STRUCTURAL_VARIANT_CALLING -> ${it}" }
        input_generation_results.bam_indexes.view { "DEBUG BAI FILES BEFORE STRUCTURAL_VARIANT_CALLING -> ${it}" }
        seqretrieval_results.genome_file.view { "DEBUG GENOME FILE BEFORE STRUCTURAL_VARIANT_CALLING -> ${it}" }

        //
        // STEP 4: Run structural variant calling
        //
        structural_results = STRUCTURAL_VARIANT_CALLING(
            input_generation_results.bam_files,
            input_generation_results.bam_indexes,
            input_generation_results.bgzip_fasta_file,
            input_generation_results.fasta_index,
            input_generation_results.bwa_index,
            input_generation_results.fasta_index_gz,
            input_generation_results.fasta_file
        )

        //
        // STEP 5: Generate final reports (e.g., MultiQC, summary statistics)
        //
        report_generation = REPORT_GENERATION()

    emit:
        multiqc_report              = data_analysis_results.multiqc_report
        structural_variants         = structural_results.delly_variants
        structural_variants_index   = structural_results.delly_variants_index
        final_reports               = report_generation.output_files
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
        ch_taxonomy_id      = Channel.value(params.taxonomy_id)
        ch_outdir           = Channel.value(params.outdir)
        ch_sequence_dir     = Channel.value(params.sequence_dir)
        ch_genome_file      = Channel.value(params.genome_file)
        ch_sequence_dir_abs = params.sequence_dir ? Channel.value(file(params.sequence_dir).toAbsolutePath()) : Channel.value('')

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
            ch_sequence_dir_abs,
            ch_sequence_dir,
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
            NFCORE_EUKAVARIZER.out.multiqc_report,
            NFCORE_EUKAVARIZER.out.structural_variants,
            NFCORE_EUKAVARIZER.out.structural_variants_index
        )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
