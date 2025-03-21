process VARIFY {
    tag "$taxonomy_id"
    conda "${moduleDir}/environment.yml"

    // Settings
    publishDir "$outdir", mode: 'copy'

    input:
        path merged_vcf
        path other_vcfs
        path survivor_vcf
        path survivor_stats
        val taxonomy_id
        val outdir

    output:
        path "index.html", emit: html_index
        path "merged_report.html", emit: html_merged
        path "survivor_report.html", emit: html_survivor

    script:
    """
    varify --merged_vcf ${merged_vcf} \\
                        \$(for vcf in ${other_vcfs}; do echo "--other_vcfs \$vcf"; done) \\
                        --survivor_vcf ${survivor_vcf} \\
                        --survivor_stats ${survivor_stats} \\
                        --output_dir .
    """
}
