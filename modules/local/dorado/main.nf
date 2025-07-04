process DORADO {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "docker.io/nanoporetech/dorado:sha268dcb4cd02093e75cdc58821f8b93719c4255ed"

    input:
    tuple val(meta), path(fast5s)

    output:
    tuple val(meta), path("${meta.id}.fastq.gz"), emit: fastq

    script:
    """
    dorado basecaller \\
        ${task.ext.model} \\
        ${task.ext.args} \\
        $fast5s \\
        | gzip > ${meta.id}.fastq.gz
    """
}
