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

include { VARIFY            } from '../../../modules/local/varify/main'

workflow REPORT_GENERATION {

    take:
        ch_taxonomy_id              // Sample taxonomy ID
        ch_outdir                   // Output directory for the reports
        vcf_list                    // List of individual VCF files
        survivor_vcf                // Merged VCF file from SURVIVOR
        survivor_stats              // SURVIVOR stats table
        bcfmerge_vcf                // Final merged VCF file from BCFtools

    main:
        // Run VARIFY to generate reports from the VCF files
        VARIFY(
            bcfmerge_vcf.map { it[1] },                              // Merged VCF file
            vcf_list.filter { it }.map { it[1] }.toList().flatten(), // Individual VCF files
            survivor_vcf.map { it[1] },                              // SURVIVOR merged VCF file
            survivor_stats.map { it[1] },                            // SURVIVOR stats table
            ch_taxonomy_id,
            ch_outdir                                                 // Output directory
        )

    emit:
        html_index      = VARIFY.out.html_index      // Summary report
        html_merged     = VARIFY.out.html_merged     // Merged VCF report
        html_survivor   = VARIFY.out.html_survivor   // SURVIVOR analysis report
}
