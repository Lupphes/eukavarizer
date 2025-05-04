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
    TIDDIT, Sniffles, CuteSV, SVABA), merges the results, and generates a structured report.
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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NFCORE_EUKAVARIZER {

    take:
        taxonomy_id
        outdir
        samplesheet
        reference_genome

    main:
        report = Channel.empty()
        sv_callers_enabled =
            params.gridss_flag || params.delly_flag || params.manta_flag ||
            params.sniffles_flag || params.cutesv_flag || params.tiddit_flag ||
            params.dysgu_flag || (params.svaba_flag && !params.bwamem2)

        log.info "🧬 SV callers enabled: ${sv_callers_enabled}"
        if (params.svaba_flag && params.bwamem2) {
            log.warn "⚠️ SVABA is not compatible with BWA-MEM2. Please use BWA-MEM v1 to run SVABA."
        }

        REFERENCE_RETRIEVAL(
            taxonomy_id,
            outdir,
            reference_genome
        )

        SEQUENCE_PROCESSOR(
            samplesheet,
            REFERENCE_RETRIEVAL.out.reference_genome_unzipped,
            REFERENCE_RETRIEVAL.out.reference_genome_bwa_index
        )

        if (sv_callers_enabled) {

            EUKAVARIZER(
                SEQUENCE_PROCESSOR.out.bam_bai,
                REFERENCE_RETRIEVAL.out.reference_genome_bgzipped,
                REFERENCE_RETRIEVAL.out.reference_genome_faidx,
                REFERENCE_RETRIEVAL.out.reference_genome_bwa_index,
                REFERENCE_RETRIEVAL.out.reference_genome_unzipped,
                REFERENCE_RETRIEVAL.out.reference_genome_minimap_index
            )

            SV_UNIFICATION (
                EUKAVARIZER.out.vcf_list,
                EUKAVARIZER.out.vcfgz_list,
                EUKAVARIZER.out.tbi_list,
                REFERENCE_RETRIEVAL.out.reference_genome_bgzipped
            )

            report = REPORT_GENERATION(
                taxonomy_id,
                outdir,
                SV_UNIFICATION.out.survivor_vcf,
                SV_UNIFICATION.out.survivor_stats,
                SV_UNIFICATION.out.bcfmerge_vcf,
                SV_UNIFICATION.out.bcfmerge_tbi,
                SV_UNIFICATION.out.bcfmerge_stats,
                EUKAVARIZER.out.vcf_list,
                EUKAVARIZER.out.tbi_list,
                REFERENCE_RETRIEVAL.out.reference_genome_unzipped
            ).report_file

        }
        else {
            log.warn "⚠️ No SV callers enabled. Skipping variant calling and report generation."
        }

    log.info "📂 Results will be saved in: ${params.outdir}"

    emit:
        multiqc_report      = Channel.empty()
        report_file         = sv_callers_enabled ? report : Channel.empty()
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
        reference_genome        = params.reference_genome ? Channel.fromPath(params.reference_genome, type: 'file', checkIfExists: true)  : []

        //
        // SUBWORKFLOW: Run initialisation tasks
        //
        PIPELINE_INITIALISATION (
            params.version,
            params.validate_params,
            params.monochrome_logs,
            args,
            params.outdir,
            params.input
        )

        //
        // WORKFLOW: Run main workflow
        //
        NFCORE_EUKAVARIZER (
            params.taxonomy_id,
            params.outdir,
            PIPELINE_INITIALISATION.out.samplesheet,
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
