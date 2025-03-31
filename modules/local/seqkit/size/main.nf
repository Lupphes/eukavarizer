process SEQKIT_SIZE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0' :
        'biocontainers/seqkit:2.9.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path(reads), path("median_size.txt"), emit: median_bp
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Compute stats once and reuse the file
    seqkit stats --tabular --all ${reads.collect { "\"${it}\"" }.join(' ')} > stats.txt

    # Extract number of sequences
    seqs=\$(awk 'NR==2 {print \$4}' stats.txt)
    if [ "\$seqs" -gt 0 ]; then
        awk 'NR==2 {printf "%d\\n", \$9}' stats.txt > median_size.txt
    else
        echo "0" > median_size.txt
    fi

    # Save version info
    echo "seqkit: \$(seqkit version 2>&1)" > versions.yml
    """
}
