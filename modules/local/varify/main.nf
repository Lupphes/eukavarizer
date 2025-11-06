process VARIFY {
    tag "$taxonomy_id"
    label 'process_low'
    cache false

    conda "${moduleDir}/environment.yml"
    container "docker.io/luppo/varify:ed5c06d6409879bfbe9a4aeeb5775409d119ba69"

    input:
    val taxonomy_id
    val outdir
    tuple val(meta1), path(survivor_vcf)
    tuple val(meta2), path(survivor_stats)
    tuple val(meta3), path(bcfmerge_vcf)
    tuple val(meta4), path(bcfmerge_tbi)
    tuple val(meta5), path(bcfmerge_stats)
    tuple val(meta6), path(reference_genome_faidx)
    tuple val(meta7), path(reference_genome_fasta)
    val(profile)

    output:
    path "*.html"             , emit: report_file
    path "genome_files/*"     , emit: genome_files
    path "versions.yml"       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${taxonomy_id}"
    """
    varify \\
        --output-dir . \\
        --bcf-vcf-file ${bcfmerge_vcf} \\
        --survivor-vcf-file ${survivor_vcf} \\
        --fasta-file ${reference_genome_fasta} \\
        --bcf-stats-file ${bcfmerge_stats} \\
        --survivor-stats-file ${survivor_stats} \\
        --report-file "${prefix}.html" \\
        --profile ${profile} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        varify: \$( varify --version 2>&1 | sed 's/varify //g' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${taxonomy_id}"
    """
    mkdir -p genome_files

    touch ${prefix}.html

    touch genome_files/${reference_genome_fasta}
    touch genome_files/${reference_genome_fasta}.fai
    touch genome_files/${bcfmerge_stats}
    touch genome_files/${survivor_stats}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        varify: \$( python -c "import varify; print(varify.__version__)" )
    END_VERSIONS
    """
}
