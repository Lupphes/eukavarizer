process SVABA_ANNOTATE {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container 'quay.io/quay.io/biocontainers/python:3.10'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.id}.reclassified.vcf"), emit: vcf

    script:
    """
    gunzip -c $vcf > input.vcf
    svaba_annotate.py input.vcf > ${meta.id}.reclassified.vcf
    """
}
