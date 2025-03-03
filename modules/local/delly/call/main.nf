process DELLY_CALL {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/delly:1.3.3--h4d20210_0"

    input:
    tuple val(meta), path(input), path(input_index)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*.{bcf,vcf.gz}"), emit: bcf
    tuple val(meta), path("*.{csi,tbi}"), emit: csi
    path "versions.yml", emit: versions

    script:
    """
    delly call \\
        --genome ${fasta} \\
        --outfile ${meta.id}.vcf.gz \\
        ${input}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        delly: \$(delly --version 2>&1 | awk '{print \$NF}')
    END_VERSIONS
    """
}
