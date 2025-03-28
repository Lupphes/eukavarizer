/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    REPORT GENERATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    This workflow handles the final report generation for the pipeline. It takes the
    merged and individual VCF files, along with the SURVIVOR analysis results, and
    generates structured HTML reports.

    The VARIFY module is used to process the VCF files and generate three reports:
    - **html_index**: Summary of the results and key metrics.
    - **html_merged**: Detailed report of the merged structural variant calls.
    - **html_survivor**: Summary of the SURVIVOR-based variant analysis.

    The reports provide a comprehensive overview of the structural variant calling results,
    allowing for easy interpretation and validation of the analysis.
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
        bcfmerge_stats
        vcf_list
        tbi_list
        reference_genome

    main:

        vcf_list = vcf_list.filter { it }
            .map { it[1] }
            .toList()

        tbi_list = tbi_list.filter { it }
            .map { it[1] }
            .toList()

        varify_meta = Channel.value([id: "varify_merge"])
        varify_input = varify_meta.combine(tbi_list.toList())

        VARIFY(
            taxonomy_id,
            outdir,
            survivor_vcf,
            survivor_stats,
            bcfmerge_vcf,
            bcfmerge_stats,
            varify_input,
            reference_genome
        )

    emit:
        report_file = VARIFY.out.report_file
}
