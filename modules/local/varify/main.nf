process VARIFY {
    tag "$taxonomy_id"
    conda "${moduleDir}/environment.yml"

    input:
        val taxonomy_id
        val outdir
        tuple val(meta1), path (survivor_vcf)
        tuple val(meta2), path (survivor_stats)
        tuple val(meta3), path (bcfmerge_vcf)
        tuple val(meta5), path (bcfmerge_stats)
        tuple val(meta6), path (varify_input)
        tuple val(meta7), path (reference_genome_bgzipped_index)

    output:
        path "report.html", emit: report_file
        path "*.png", emit: report_images

    script:
    """
    varify \\
        --output-dir . \\
        --bfc-vcf-file ${bcfmerge_vcf} \\
        --survivor-vcf-file ${survivor_vcf} \\
        --fasta-file ${reference_genome_bgzipped_index} \\
        --bfc-stats-file ${bcfmerge_stats} \\
        --survivor-stats-file ${survivor_stats} \\
        --report-file "report.html"

    """
}
// \$(for vcf in ${vcf_list}; do echo "--sample-vcf-files \$vcf"; done) \\
