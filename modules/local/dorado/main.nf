process DORADO {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "docker.io/nanoporetech/dorado:sha268dcb4cd02093e75cdc58821f8b93719c4255ed"

    input:
    tuple val(meta), path(fast5s)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def model = task.ext.model ?: 'dna_r9.4.1_e8_hac@v3.3'
    """
    dorado basecaller \\
        ${model} \\
        ${args} \\
        ${fast5s} \\
        | gzip > ${prefix}.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: \$( dorado --version 2>&1 | sed 's/dorado //g' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf "@read1\\nACGTACGT\\n+\\n!!!!!!!!\\n" | gzip > ${prefix}.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: \$( dorado --version 2>&1 | sed 's/dorado //g' )
    END_VERSIONS
    """
}
