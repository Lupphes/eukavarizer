process GRIDSS_ANNOTATE {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"

    container 'quay.io/biocontainers/bioconductor-structuralvariantannotation:1.22.0--r44hdfd78af_0'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.id}.reclassified.vcf"), emit: vcf
    path("${meta.id}.simple.bed"), emit: bed

    script:
    """
    gunzip -c $vcf > input.vcf
    simple-event-annotation.R input.vcf ${meta.id}.reclassified.vcf ${meta.id}.simple.bed
    """
}
