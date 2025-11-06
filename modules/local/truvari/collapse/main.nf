process TRUVARI_COLLAPSE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/truvari:5.3.0--pyhdfd78af_0':
        'biocontainers/truvari:5.3.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*.collapsed.vcf.gz")     , emit: collapsed_vcf
    tuple val(meta), path("*.removed.vcf")          , emit: removed_vcf, optional: true
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def keep_removed = task.ext.keep_removed ? "--removed-output ${prefix}.removed.vcf" : ""

    """
    truvari collapse \\
        --input ${vcf} \\
        --output ${prefix}.collapsed.vcf.gz \\
        --reference ${fasta} \\
        ${keep_removed} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        truvari: \$(echo \$(truvari version 2>&1) | sed 's/^Truvari v//' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip > ${prefix}.collapsed.vcf.gz
    touch ${prefix}.removed.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        truvari: \$(echo \$(truvari version 2>&1) | sed 's/^Truvari v//' )
    END_VERSIONS
    """
}