process BIODBCORE_REFSEQ {
    tag "$taxonomy_id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "docker.io/luppo/biodbcore:1b8dd90e307a6937b09b9d7088f3b8ec5e8dacc8"

    input:
    val  taxonomy_id
    val  outdir
    path reference_genome

    output:
    path "*_refseq_results.json"         , emit: json
    path "*.{fna.gz,fa.gz}"              , emit: reference_genome, optional: true
    path "assembly_summary_refseq.parquet", emit: refseq_parquet, optional: true
    path "assembly_summary_refseq.txt"   , emit: refseq_summary, optional: true
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${taxonomy_id}"
    def reference_arg = reference_genome ? "--reference_genome ${reference_genome}" : ''
    """
    if [ -n "${reference_genome}" ]; then
        echo "Reference genome already exists at: ${reference_genome} for taxonomy ID: ${taxonomy_id}"
    else
        echo "Downloading reference genome for taxonomy ID: ${taxonomy_id}"
    fi

    biodbcore \\
        --mode refseq \\
        --taxonomy_id ${taxonomy_id} \\
        --outdir . \\
        ${reference_arg} \\
        ${args}

    # Rename output with configurable prefix
    mv refseq_results.json ${prefix}_refseq_results.json

    cat ${prefix}_refseq_results.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biodbcore: \$( biodbcore --version 2>&1 | sed 's/biodbcore //g' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${taxonomy_id}"
    """
    touch ${prefix}_refseq_results.json

    # Create mock reference genome if not provided
    if [ -z "${reference_genome}" ]; then
        echo ">mock_chromosome_1" | gzip > ${prefix}_mock_genome.fna.gz
    fi

    # Create mock assembly files
    touch assembly_summary_refseq.txt
    touch assembly_summary_refseq.parquet

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biodbcore: \$( biodbcore --version 2>&1 | sed 's/biodbcore //g' )
    END_VERSIONS
    """
}
