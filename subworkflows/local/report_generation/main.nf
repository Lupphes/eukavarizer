/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    REPORT GENERATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    This workflow handles the final report generation for the pipeline. It takes the
    merged and individual VCF files, along with the SURVIVOR analysis results, and
    generates structured HTML reports.

    The VARIFY module is used to process the VCF files and generate three reports:
    - **html_index** – Summary of the results and key metrics.
    - **html_merged** – Detailed report of the merged structural variant calls.
    - **html_survivor** – Summary of the SURVIVOR-based variant analysis.

    Outputs:
    - `report_file`     – Main report (`report.html`).
    - `report_images`   – PNG images (figures and plots).
    - `report_html`     – Additional HTML report files from the plots directory.
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

    main:
        vcf_list_cleaned = vcf_list
            .filter { it != null }
            .map { it[1] }

        tbi_list_cleaned = tbi_list
            .filter { it != null }
            .map { it[1] }

        varify_meta = Channel.value([id: "varify_merge"])
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
            workflow.profile
        )

    emit:
        report_file = VARIFY.out.report_file
}
