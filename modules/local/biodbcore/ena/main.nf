process BIODBCORE_ENA {
    tag "$taxonomy_id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "luppo/biodbcore:latest"

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
        path "ena_results.json", emit: ena_results
        path "samplesheet.csv", emit: samplesheet

    script:
    """
    echo "No local sequencing data found. Downloading sequencing reads for taxonomy ID: $taxonomy_id"

    biodbcore --mode ena \
        --taxonomy_id $taxonomy_id \
        --library_strategy $library_strategy \
        --instrument_platform $instrument_platform \
        --minimum_coverage $minimum_coverage \
        --maximum_coverage $maximum_coverage \
        --max_results $max_results \
        \$( [ -n "$assembly_quality" ] && echo "--assembly_quality $assembly_quality" ) \
        --outdir . \
        --genome_size_ungapped $genome_size_ungapped


    echo    "Sequencing data retrieval complete."
    cat     "ena_results.json"
    """
}
