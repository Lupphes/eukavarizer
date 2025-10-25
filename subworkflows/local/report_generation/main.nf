/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: REPORT_GENERATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Generates comprehensive HTML reports for structural variant analysis results
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Description:
        This subworkflow creates detailed HTML reports summarizing structural variant calls,
        quality metrics, and analysis statistics. It processes merged and individual VCF files
        to generate visualizations and summaries using the VARIFY reporting tool.

    Processing Steps:
        1. Collect and clean VCF and index files from all callers
        2. Collect alignment statistics from Samtools
        3. Generate comprehensive report with VARIFY module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Input Channels (take):
        taxonomy_id               val(id)                NCBI taxonomy identifier
        outdir                    path(dir)              Output directory path
        survivor_vcf              tuple(meta, vcf)       SURVIVOR merged VCF
        survivor_stats            tuple(meta, stats)     SURVIVOR statistics
        bcfmerge_vcf              tuple(meta, vcf.gz)    BCFtools merged VCF
        bcfmerge_tbi              tuple(meta, tbi)       BCFtools VCF index
        bcfmerge_stats            tuple(meta, stats)     BCFtools statistics
        vcf_list                  tuple(meta, vcf)       Individual caller VCFs
        tbi_list                  tuple(meta, tbi)       Individual caller indexes
        reference_genome          tuple(meta, fasta)     Reference genome
        samtools_stats            tuple(meta, stats)     Alignment statistics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Output Channels (emit):
        report_file               path(html)             Main analysis report
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Author:   Ondřej Sloup (Lupphes)
    Contact:  ondrej.sloup@protonmail.com
    GitHub:   @Lupphes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { VARIFY } from '../../../modules/local/varify/main'

workflow REPORT_GENERATION {

    take:
        taxonomy_id
        outdir
        survivor_vcf
        survivor_stats
        bcfmerge_vcf
        bcfmerge_tbi
        bcfmerge_stats
        vcf_list
        tbi_list
        reference_genome
        samtools_stats

    main:
        ch_versions = channel.empty()

        vcf_list_cleaned = vcf_list
            .filter { tuple -> tuple != null }
            .map { tuple -> tuple[1] }

        _tbi_list_cleaned = tbi_list
            .filter { tuple -> tuple != null }
            .map { tuple -> tuple[1] }

        samtools_stats_cleaned = samtools_stats
            .map { tuple -> tuple[1] }
            .collect()

        varify_meta = channel.value([id: "varify_merge"])
        varify_input = varify_meta.combine(vcf_list_cleaned)

        VARIFY(
            taxonomy_id,
            outdir,
            survivor_vcf,
            survivor_stats,
            bcfmerge_vcf,
            bcfmerge_tbi,
            bcfmerge_stats,
            varify_input,
            reference_genome,
            workflow.profile,
            samtools_stats_cleaned
        )

        ch_versions = ch_versions.mix(VARIFY.out.versions)

    emit:
        report_file = VARIFY.out.report_file
        versions = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
