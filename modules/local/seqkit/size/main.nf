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
    # Extract median size using seqkit and awk
    seqkit seq -s ${reads} | awk '{print length}' | sort -n | awk '{a[NR]=\$1} END {if (NR%2==1) {print a[(NR+1)/2]} else {print (a[NR/2] + a[NR/2+1])/2}}' > median_size.txt

    # Record seqkit version info
    echo "seqkit: \$(seqkit version 2>&1)" > versions.yml
    """
}
