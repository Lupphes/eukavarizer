process GZIP {
    tag "${input_file}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data'
        : 'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8'}"

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("${gzip_file}"), emit: gzip
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def name = input_file.getBaseName()
    gzip_file = "${name}.gz"
    """
    gzip \\
        -c \\
        ${args} \\
        ${input_file} \\
        > ${gzip_file}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gzip: \$(echo \$(gzip --version 2>&1) | sed 's/^.*(gzip) //; s/ Copyright.*\$//')
    END_VERSIONS
    """

    stub:
    def name = input_file.getBaseName()
    gzip_file = "${name}.gz"
    """
    touch ${gzip_file}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gzip: \$(echo \$(gzip --version 2>&1) | sed 's/^.*(gzip) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}
