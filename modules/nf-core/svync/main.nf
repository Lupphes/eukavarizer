process SVYNC {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/svync:0.3.0--h9ee0642_0':
        'biocontainers/svync:0.3.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi), path(config)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf, optional: true
    tuple val(meta), path("*.tbi")   , emit: tbi, optional: true
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def args2  = task.ext.args2 ?: ''
    def args3  = task.ext.args3 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if ("$vcf" == "${prefix}.vcf.gz") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    """
    # Function to read VCF (handles both compressed and uncompressed)
    read_vcf() {
        if [[ "$vcf" == *.gz ]]; then
            zcat "\$1"
        else
            cat "\$1"
        fi
    }

    # Check if VCF has samples by counting columns in #CHROM header
    # Sites-only VCFs have 8 columns (no FORMAT/sample), sample VCFs have 9+ columns
    HEADER_LINE=\$(read_vcf $vcf | grep "^#CHROM" | head -n 1)
    COLUMN_COUNT=\$(echo "\$HEADER_LINE" | awk '{print NF}')

    # Check if VCF has any variants (non-header lines)
    VARIANT_COUNT=\$(read_vcf $vcf | grep -v "^#" | wc -l || echo "0")

    if [[ "\$COLUMN_COUNT" -le 8 ]]; then
        # Sites-only VCF (no sample columns) - svync cannot process these
        # Do not create output files - flow will stop here (like SAMPLE_REHEADER does)
        echo "WARNING: Input VCF is sites-only (no sample columns). Skipping svync processing." >&2

    elif [[ "\$VARIANT_COUNT" -eq 0 ]]; then
        # Empty VCF (no variants) - svync cannot process
        # Do not create output files - flow will stop here
        echo "WARNING: Input VCF has no variants. Skipping svync processing." >&2

    else
        # Normal VCF with samples and variants - process with svync
        svync \\
            $args \\
            --config $config \\
            --input $vcf \\
            | bgzip --threads $task.cpus $args2 > ${prefix}.vcf.gz \\
            && tabix $args3 ${prefix}.vcf.gz
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svync: \$(svync --version | sed 's/svync version //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    if ("$vcf" == "${prefix}.vcf.gz") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    """
    echo | gzip -n > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svync: \$(svync --version | sed 's/svync version //')
    END_VERSIONS
    """
}
