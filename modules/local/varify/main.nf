process VARIFY {
    tag "$taxonomy_id"
    label 'process_low'
    cache false

    conda "${moduleDir}/environment.yml"
    container "docker.io/luppo/varify:ac8b977976d91caf2495d1e02cb135073fac722b"

    input:
    val taxonomy_id
    val outdir
    tuple val(meta1), path(survivor_vcf)
    tuple val(meta2), path(survivor_stats)
    tuple val(meta3), path(bcfmerge_vcf)
    tuple val(meta4), path(bcfmerge_tbi)
    tuple val(meta5), path(bcfmerge_stats)
    tuple val(meta6), path(vcf_list_cleaned)
    tuple val(meta7), path(reference_genome_bgzipped_index)
    val(profile)
    path(samtools_stats)

    output:
    path "*.html"       , emit: report_file
    path "plots/*.png"  , emit: report_images
    path "plots/*.html" , emit: report_html
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${taxonomy_id}"
    def samtools_stats_files = samtools_stats ? samtools_stats.collect { "--samtools-stats-file ${it}" }.join(' ') : ''
    """
    varify \\
        --output-dir . \\
        --bcf-vcf-file ${bcfmerge_vcf} \\
        --survivor-vcf-file ${survivor_vcf} \\
        --fasta-file ${reference_genome_bgzipped_index} \\
        --bcf-stats-file ${bcfmerge_stats} \\
        --survivor-stats-file ${survivor_stats} \\
        --report-file "${prefix}.html" \\
        --profile ${profile} \\
        # TODO: Add samtools parameter
        # ${samtools_stats_files} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        varify: \$( varify --version 2>&1 | sed 's/varify //g' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${taxonomy_id}"
    """
    mkdir -p plots
    touch ${prefix}.html
    touch plots/dummy.png
    touch plots/dummy.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        varify: \$( varify --version 2>&1 | sed 's/varify //g' )
    END_VERSIONS
    """
}
