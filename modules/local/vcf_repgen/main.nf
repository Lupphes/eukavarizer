process VCF_REPGEN {
    tag "$taxonomy_id"

    conda "${moduleDir}/environment.yml"

    publishDir "${outdir}/reports", mode: 'copy'

    input:
        path merged_vcf
        path other_vcfs
        path survivor_vcf
        path survivor_stats
        val taxonomy_id
        val outdir

    output:
        path "${taxonomy_id}/report.html", emit: html

    script:
    """
    vcf_gen.py --merged_vcf ${merged_vcf} \\
                        \$(for vcf in ${other_vcfs}; do echo "--other_vcfs \$vcf"; done) \\
                        --survivor_vcf ${survivor_vcf} \\
                        --survivor_stats ${survivor_stats} \\
                        --output_dir ${outdir}/${taxonomy_id}
    """
}
