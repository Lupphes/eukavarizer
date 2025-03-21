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

nextflow.enable.dsl = 2
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { EUKAVARIZER               } from './workflows/eukavarizer'

include { REFERENCE_RETRIEVAL       } from './subworkflows/local/reference_retrieval'
include { SEQUENCE_PROCESSOR        } from './subworkflows/local/sequence_processor'
include { SV_UNIFICATION            } from './subworkflows/local/sv_unification'
include { REPORT_GENERATION         } from './subworkflows/local/report_generation'

// include { SEQRETRIEVAL                  } from './workflows/seqretrieval'
// include { ANALYSIS_FAST_MULTI           } from './workflows/analysis_fast_multi'
// include { INPUT_GENERATION              } from './workflows/input_generation'
// include { REPORT_GENERATION  as OUTDATED            } from './workflows/report_generation'
include { PIPELINE_INITIALISATION   } from './subworkflows/local/utils_nfcore_eukavarizer_pipeline'
include { PIPELINE_COMPLETION       } from './subworkflows/local/utils_nfcore_eukavarizer_pipeline'
// include { STRUCTURAL_VARIANT_CALLING    } from './workflows/structural_variant_calling'

params.gridss_flag                  = true // Tech Debt: Gridss is always enabled

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NFCORE_EUKAVARIZER {

    take:
        ch_taxonomy_id
        ch_outdir
        ch_sequences_abs_dir
        ch_sequence_dir
        ch_reference_genome
        ch_max_distance_breakpoints
        ch_min_supporting_callers
        ch_account_for_type
        ch_account_for_sv_strands
        ch_estimate_distanced_by_sv_size
        ch_min_sv_size

    main:

        REFERENCE_RETRIEVAL(
            ch_taxonomy_id,
            ch_outdir,
            ch_reference_genome
        )

        SEQUENCE_PROCESSOR(
            ch_taxonomy_id,
            ch_outdir,
            ch_sequences_abs_dir,
            REFERENCE_RETRIEVAL.out.reference_genome_ungapped_size,
            REFERENCE_RETRIEVAL.out.reference_genome_unzipped,
            REFERENCE_RETRIEVAL.out.reference_genome_bgzipped,
            REFERENCE_RETRIEVAL.out.reference_genome_bwa_index
        )

        EUKAVARIZER(
            SEQUENCE_PROCESSOR.out.fastq_bam,
            SEQUENCE_PROCESSOR.out.fastq_bam_indexes,
            REFERENCE_RETRIEVAL.out.reference_genome_bgzipped,
            REFERENCE_RETRIEVAL.out.reference_genome_faidx,
            REFERENCE_RETRIEVAL.out.reference_genome_bwa_index,
            REFERENCE_RETRIEVAL.out.reference_genome_bgzipped_faidx,
            REFERENCE_RETRIEVAL.out.reference_genome_unzipped
        )

        SV_UNIFICATION (
        )

        // REPORT_GENERATION()

        // //
        // // STEP 1: Retrieve Sequences & Reference Genome
        // //
        // SEQRETRIEVAL(
        //     ch_taxonomy_id,
        //     ch_outdir,
        //     ch_sequences_abs_dir,
        //     ch_reference_genome,
        //     debug_flag
        // )

        //
        // STEP 2: Run Main Analysis Pipeline
        //
        // data_analysis_results = ANALYSIS_FAST_MULTI(
        //     seqretrieval_results.grouped_fastqs,
        //     debug_flag,
        //     fastqc_flag,
        //     multiqc_flag
        // )

        // //
        // // STEP 3: Preprocess and generate INPUT for Structural Analysis
        // // with alignments
        // //
        // view("🔹 Checking inputs for INPUT_GENERATION")

        // input_generation_results = INPUT_GENERATION(
        //     seqretrieval_results.genome_file,
        //     data_analysis_results.processed_files,
        //     debug_flag
        // )

        // //
        // // STEP 4: Run structural variant calling
        // //
        // structural_results = STRUCTURAL_VARIANT_CALLING(
        //     input_generation_results.bam_files,
        //     input_generation_results.bam_indexes,
        //     input_generation_results.bgzip_fasta_file,
        //     input_generation_results.fasta_index,
        //     input_generation_results.bwa_index,
        //     input_generation_results.fasta_index_gz,
        //     input_generation_results.fasta_file,
        //     delly_flag,
        //     manta_flag,
        //     gridss_flag,
        //     dysgu_flag,
        //     tiddit_flag,
        //     svaba_flag,
        //     sniffles_flag,
        //     cutesv_flag,
        //     debug_flag
        // )

        // //
        // // STEP 5: Generate final reports (e.g., MultiQC, summary statistics)
        // //
        // report_generation = REPORT_GENERATION(
        //     seqretrieval_results.genome_file,
        //     structural_results.delly_variants,
        //     structural_results.delly_variants_bgzipped,
        //     structural_results.delly_variants_index,
        //     structural_results.delly_variants_index_bgzipped,
        //     structural_results.manta_small_variants,
        //     structural_results.manta_small_variants_index,
        //     structural_results.manta_small_variants_bgzipped,
        //     structural_results.manta_candidate_variants,
        //     structural_results.manta_candidate_variants_index,
        //     structural_results.manta_candidate_variants_bgzipped,
        //     structural_results.manta_diploid_variants,
        //     structural_results.manta_diploid_variants_index,
        //     structural_results.manta_diploid_variants_bgzipped,
        //     structural_results.gridss_variants,
        //     structural_results.gridss_variants_index,
        //     structural_results.gridss_variants_bgzipped,
        //     structural_results.dysgu_variants,
        //     structural_results.dysgu_variants_index,
        //     structural_results.dysgu_variants_bgzipped,
        //     structural_results.tiddit_variants,
        //     structural_results.tiddit_ploidy,
        //     structural_results.tiddit_variants_index,
        //     structural_results.tiddit_variants_bgzipped,
        //     structural_results.cutesv_variants,
        //     structural_results.cutesv_variants_index,
        //     structural_results.cutesv_variants_bgzipped,
        //     delly_flag,
        //     manta_flag,
        //     gridss_flag,
        //     dysgu_flag,
        //     tiddit_flag,
        //     svaba_flag,
        //     sniffles_flag,
        //     cutesv_flag,
        //     debug_flag,
        //     ch_max_distance_breakpoints,
        //     ch_min_supporting_callers,
        //     ch_account_for_type,
        //     ch_account_for_sv_strands,
        //     ch_estimate_distanced_by_sv_size,
        //     ch_min_sv_size,
        //     ch_taxonomy_id,
        //     ch_outdir
        // )

    emit:
        multiqc_report              = "data_analysis_results.multiqc_report"
        structural_variants         = "structural_results.delly_variants"
        structural_variants_index   = "structural_results.delly_variants_index"
        html_index                  = "report_generation.html_index"
        html_merged                 = "report_generation.html_merged"
        html_survivor               = "report_generation.html_survivor"
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
        //
        // Input channels for the workflow
        //
        ch_taxonomy_id              = Channel.value(params.taxonomy_id)
        ch_outdir                   = Channel.value(params.outdir)
        ch_sequence_dir             = Channel.value(params.sequence_dir)
        ch_reference_genome         = Channel.value(params.genome_file)
        ch_sequence_dir_abs = params.sequence_dir ? Channel.value(file(params.sequence_dir).toAbsolutePath()) : Channel.value('')



        // Report Options
        ch_max_distance_breakpoints         = params.max_distance_breakpoints ? Channel.value(params.max_distance_breakpoints) : Channel.value(1000)
        ch_min_supporting_callers           = params.min_supporting_callers ? Channel.value(params.min_supporting_callers) : Channel.value(1)
        ch_account_for_type                 = params.account_for_type ? Channel.value(params.account_for_type) : Channel.value(1)
        ch_account_for_sv_strands           = params.account_for_sv_strands ? Channel.value(params.account_for_sv_strands) : Channel.value(0)
        ch_estimate_distanced_by_sv_size    = params.estimate_distanced_by_sv_size ? Channel.value(params.estimate_distanced_by_sv_size) : Channel.value(0)
        ch_min_sv_size                      = params.min_sv_size ? Channel.value(params.min_sv_size) : Channel.value(30)

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
            ch_reference_genome,
            ch_max_distance_breakpoints,
            ch_min_supporting_callers,
            ch_account_for_type,
            ch_account_for_sv_strands,
            ch_estimate_distanced_by_sv_size,
            ch_min_sv_size
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
            NFCORE_EUKAVARIZER.out.structural_variants_index,
            NFCORE_EUKAVARIZER.out.html_index,
            NFCORE_EUKAVARIZER.out.html_merged,
            NFCORE_EUKAVARIZER.out.html_survivor
        )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
