process GRIDSS_ANNOTATE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'quay.io/biocontainers/bioconductor-structuralvariantannotation:1.22.0--r44hdfd78af_0'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.reclassified.vcf"), emit: vcf
    path "*.simple.bed"                        , emit: bed
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gunzip -c ${vcf} > input.vcf

    # Check if VCF has any variants (non-header lines)
    # Use || true to prevent grep from causing pipeline failure when no matches found
    VARIANT_COUNT=\$(grep -v "^#" input.vcf | wc -l || true)

    if [ "\${VARIANT_COUNT}" -eq 0 ]; then
        # No variants - create empty output files with proper headers
        grep "^#" input.vcf > ${prefix}.reclassified.vcf
        touch ${prefix}.simple.bed
        echo "No variants found in VCF. Created empty output files." >&2
    else
        # Process normally with R script
        simple-event-annotation.R ${args} input.vcf ${prefix}.reclassified.vcf ${prefix}.simple.bed
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$( R --version 2>&1 | head -n1 | sed 's/R version //; s/ .*//' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.reclassified.vcf
    touch ${prefix}.simple.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$( R --version 2>&1 | head -n1 | sed 's/R version //; s/ .*//' )
    END_VERSIONS
    """
}
