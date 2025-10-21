process BIODBCORE_ENA {
    tag "$taxonomy_id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "docker.io/luppo/biodbcore:1b8dd90e307a6937b09b9d7088f3b8ec5e8dacc8"

    input:
    val taxonomy_id
    val genome_size_ungapped
    val outdir
    val library_strategy
    val instrument_platform
    val minimum_coverage
    val maximum_coverage
    val max_results
    val assembly_quality

    output:
    path "*_ena_results.json", emit: ena_results
    path "*_samplesheet.csv"  , emit: samplesheet
    path "versions.yml"       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${taxonomy_id}"
    def assembly_quality_arg = assembly_quality ? "--assembly_quality ${assembly_quality}" : ''
    """
    echo "No local sequencing data found. Downloading sequencing reads for taxonomy ID: ${taxonomy_id}"

    biodbcore \\
        --mode ena \\
        --taxonomy_id ${taxonomy_id} \\
        --library_strategy ${library_strategy} \\
        --instrument_platform ${instrument_platform} \\
        --minimum_coverage ${minimum_coverage} \\
        --maximum_coverage ${maximum_coverage} \\
        --max_results ${max_results} \\
        ${assembly_quality_arg} \\
        --outdir . \\
        --genome_size_ungapped ${genome_size_ungapped} \\
        ${args}

    # Rename outputs with configurable prefix
    mv ena_results.json ${prefix}_ena_results.json
    mv samplesheet.csv ${prefix}_samplesheet.csv

    echo "Sequencing data retrieval complete."
    cat ${prefix}_ena_results.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biodbcore: \$( biodbcore --version 2>&1 | sed 's/biodbcore //g' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${taxonomy_id}"
    """
    touch ${prefix}_ena_results.json
    touch ${prefix}_samplesheet.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biodbcore: \$( biodbcore --version 2>&1 | sed 's/biodbcore //g' )
    END_VERSIONS
    """
}
