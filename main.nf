#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/eukavarizer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/eukavarizer
    Website: https://nf-co.re/eukavarizer
    Slack  : https://nfcore.slack.com/channels/eukavarizer
----------------------------------------------------------------------------------------
    This pipeline is designed to process and analyse structural variants (SVs) in
    eukaryotic genomes. It supports multiple SV callers (DELLY, Manta, GRIDSS, Dysgu,
    TIDDIT, Sniffles, CuteSV), merges the results, and generates a structured report.
    The pipeline follows nf-core best practices for modularity and reproducibility.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

include { EUKAVARIZER               } from './workflows/eukavarizer'

include { REFERENCE_RETRIEVAL       } from './subworkflows/local/reference_retrieval'
include { SEQUENCE_PROCESSOR        } from './subworkflows/local/sequence_processor'
include { SV_UNIFICATION            } from './subworkflows/local/sv_unification'
include { REPORT_GENERATION         } from './subworkflows/local/report_generation'

include { PIPELINE_INITIALISATION   } from './subworkflows/local/utils_nfcore_eukavarizer_pipeline'
include { PIPELINE_COMPLETION       } from './subworkflows/local/utils_nfcore_eukavarizer_pipeline'

// Sequence Processor Parameters
params.minimap2_threshold               = params.minimap2_threshold ?: 300

// Eukavarizer Pipeline Parameters
params.gridss_flag                      = true // Tech Debt: Gridss is always enabled

// Report Generation Parameters
params.max_distance_breakpoints         = params.max_distance_breakpoints ?: 1000
params.min_supporting_callers           = params.min_supporting_callers ?: 1
params.account_for_type                 = params.account_for_type ?: 1
params.account_for_sv_strands           = params.account_for_sv_strands ?: 0
params.estimate_distanced_by_sv_size    = params.estimate_distanced_by_sv_size ?: 0
params.min_sv_size                      = params.min_sv_size ?: 30

// Survivor Filter Parameters
params.min_sv_size_filter               = params.min_sv_size_filter ?: 50
params.max_sv_size_filter               = params.max_sv_size_filter ?: 100000
params.min_allele_freq_filter           = params.min_allele_freq_filter ?: 0.01
params.min_num_reads_filter             = params.min_num_reads_filter ?: 3

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
            EUKAVARIZER.out.vcf_list,
            EUKAVARIZER.out.vcfgz_list,
            EUKAVARIZER.out.tbi_list
        )

        REPORT_GENERATION(
            ch_taxonomy_id,
            ch_outdir,
            EUKAVARIZER.out.vcf_list,
            SV_UNIFICATION.out.survivor_vcf,
            SV_UNIFICATION.out.survivor_stats,
            SV_UNIFICATION.out.bcfmerge_vcf
        )

    emit:
        multiqc_report      = "data_analysis_results.multiqc_report"
        html_index          = REPORT_GENERATION.out.html_index
        html_merged         = REPORT_GENERATION.out.html_merged
        html_survivor       = REPORT_GENERATION.out.html_survivor
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
        taxonomy_id              = Channel.value(params.taxonomy_id)
        outdir                   = Channel.value(params.outdir)
        sequence_dir             = Channel.value(params.sequence_dir)
        reference_genome         = Channel.value(params.genome_file)
        sequence_dir_abs         = params.sequence_dir ? Channel.value(file(params.sequence_dir).toAbsolutePath()) : Channel.value('')

        //
        // SUBWORKFLOW: Run initialisation tasks
        //
        PIPELINE_INITIALISATION (
            params.version,
            params.validate_params,
            params.monochrome_logs,
            args,
            outdir
        )

        //
        // WORKFLOW: Run main workflow
        //
        NFCORE_EUKAVARIZER (
            taxonomy_id,
            outdir,
            sequence_dir_abs,
            sequence_dir,
            reference_genome
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
