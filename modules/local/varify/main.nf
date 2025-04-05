process VARIFY {
    tag "$taxonomy_id"
    label 'process_low'
    cache false

    conda "${moduleDir}/environment.yml"
    container "luppo/varify:latest"

    input:
        val taxonomy_id
        val outdir
        tuple val(meta1), path (survivor_vcf)
        tuple val(meta2), path (survivor_stats)
        tuple val(meta3), path (bcfmerge_vcf)
        tuple val(meta4), path (bcfmerge_tbi)
        tuple val(meta5), path (bcfmerge_stats)
        tuple val(meta6), path (vcf_list_cleaned)
        tuple val(meta7), path (reference_genome_bgzipped_index)
        val(profile)

    output:
        path "report.html", emit: report_file
        path "plots/*.png", emit: report_images
        path "plots/*.html", emit: report_html

    script:
    """
    varify \\
        --output-dir . \\
        --bcf-vcf-file ${bcfmerge_vcf} \\
        --survivor-vcf-file ${survivor_vcf} \\
        --fasta-file ${reference_genome_bgzipped_index} \\
        --bcf-stats-file ${bcfmerge_stats} \\
        --survivor-stats-file ${survivor_stats} \\
        --report-file "report.html" \\
        --profile ${profile}
    """
}
