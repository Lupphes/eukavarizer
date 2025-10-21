process SVABA_ANNOTATE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'quay.io/biocontainers/python:3.10'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.reclassified.vcf"), emit: vcf
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gunzip -c ${vcf} > input.vcf
    svaba_annotate.py ${args} input.vcf > ${prefix}.reclassified.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version 2>&1 | sed 's/Python //g' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.reclassified.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version 2>&1 | sed 's/Python //g' )
    END_VERSIONS
    """
}
