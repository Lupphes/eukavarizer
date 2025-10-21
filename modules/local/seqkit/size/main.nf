process SEQKIT_SIZE {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path(reads), path("*.txt"), emit: median_bp
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Compute quick average read length from first 1000 reads (4000 lines)
    ( zcat ${reads.collect { "\"${it}\"" }.join(' ')} || true ) | \\
        head -n 4000 | \\
        awk 'NR % 4 == 2 { total += length(\$0); count++ } END { if (count > 0) printf "%d\\n", total / count; else print 0 }' > ${prefix}_median_size.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version 2>&1 | head -n1 | sed 's/.*Awk //; s/,.*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "0" > ${prefix}_median_size.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version 2>&1 | head -n1 | sed 's/.*Awk //; s/,.*//')
    END_VERSIONS
    """
}
