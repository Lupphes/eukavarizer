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

// Default SURVIVOR filter parameters (when no profile is set)
params.sur_max_distance_breakpoints         = params.sur_max_distance_breakpoints       ?: 2000  // Allow a larger gap between breakpoints, useful for large SVs
params.sur_min_supporting_callers           = params.sur_min_supporting_callers         ?: 1     // Require at least 1 supporting caller
params.sur_account_for_type                 = params.sur_account_for_type               ?: 1     // Consider all types of SVs (deletions, duplications, etc.)
params.sur_account_for_sv_strands           = params.sur_account_for_sv_strands         ?: 0     // Don't account for strands (optional, can be adjusted based on needs)
params.sur_estimate_distanced_by_sv_size    = params.sur_estimate_distanced_by_sv_size  ?: 0     // Don't estimate distance based on SV size
params.sur_min_sv_size                      = params.sur_min_sv_size                    ?: 50    // Minimum SV size of 50 bp (small enough for small variants)

params.sur_min_sv_size_filter               = params.sur_min_sv_size_filter             ?: 100   // Filter out very small variants (set to 100 bp)
params.sur_max_sv_size_filter               = params.sur_max_sv_size_filter             ?: 50000 // Maximum SV size allowed in the final output (50 kbp as a reasonable upper limit)
params.sur_min_allele_freq_filter           = params.sur_min_allele_freq_filter         ?: 0.05  // Set to 5% allele frequency for minimum threshold
params.sur_min_num_reads_filter             = params.sur_min_num_reads_filter           ?: 3     // Minimum 3 supporting reads (default for most datasets)



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NFCORE_EUKAVARIZER {

    take:
        taxonomy_id
        outdir
        sequences_abs_dir
        reference_genome

    main:

        REFERENCE_RETRIEVAL(
            taxonomy_id,
            outdir,
            reference_genome
        )

        SEQUENCE_PROCESSOR(
            taxonomy_id,
            outdir,
            sequences_abs_dir,
            REFERENCE_RETRIEVAL.out.reference_genome_ungapped_size,
            REFERENCE_RETRIEVAL.out.reference_genome_unzipped,
            REFERENCE_RETRIEVAL.out.reference_genome_bgzipped,
            REFERENCE_RETRIEVAL.out.reference_genome_bwa_index
        )

        if (params.gridss_flag || params.delly_flag || params.manta_flag ||
            params.sniffles_flag || params.cutesv_flag || params.tiddit_flag ||
            params.dysgu_flag || params.svaba_flag) {

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
                EUKAVARIZER.out.tbi_list,
                REFERENCE_RETRIEVAL.out.reference_genome_bgzipped
            )

            REPORT_GENERATION(
                taxonomy_id,
                outdir,
                SV_UNIFICATION.out.survivor_vcf,
                SV_UNIFICATION.out.survivor_stats,
                SV_UNIFICATION.out.bcfmerge_vcf,
                SV_UNIFICATION.out.bcfmerge_stats,
                EUKAVARIZER.out.vcf_list,
                EUKAVARIZER.out.tbi_list,
                REFERENCE_RETRIEVAL.out.reference_genome_unzipped,
            )

        }
        else {
            log.warning "No SV callers enabled."
        }

    emit:
        multiqc_report      = "SEQUENCE_PROCESSOR.out.multiqc_report"
        report_file         = REPORT_GENERATION.out.report_file
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
        reference_genome         = params.reference_genome ? Channel.fromPath(params.reference_genome, type: 'file', checkIfExists: true)  : []
        sequence_dir         = params.sequence_dir ? Channel.fromPath(params.sequence_dir, type: 'dir', checkIfExists: true) : []

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
