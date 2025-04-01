process SEQKIT_SIZE {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path(reads), path("median_size.txt"), emit: median_bp
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Compute quick average read length from first 1000 reads (4000 lines)
    ( zcat ${reads.collect { "\"${it}\"" }.join(' ')} || true ) | \\
        head -n 4000 | \\
        awk 'NR % 4 == 2 { total += length(\$0); count++ } END { if (count > 0) printf "%d\\n", total / count; else print 0 }' > median_size.txt

    echo "seqkit: not used (fast method)" > versions.yml
    """
}
