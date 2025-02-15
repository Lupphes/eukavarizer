process SEQUENCE_RETRIEVER {
    // Docker
    // container 'eukalyzer/seq_retriver:latest'
    // containerOptions "-v ${params.outdir}/data:/data"

    errorStrategy 'terminate'
    publishDir "${params.outdir}/sequences", mode: 'copy'

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
        val sequences_dir

    output:
        path "$outdir/sequence_results.json", emit: sequence_json

    script:
    """
    echo "Downloading sequencing reads for taxonomy ID: $taxonomy_id"

    python /app/main.py --mode ena \
        --taxonomy_id $taxonomy_id \
        --library_strategy '${library_strategy.join(",")}' \
        --instrument_platform '${instrument_platform.join(",")}' \
        --minimum_coverage $minimum_coverage \
        --maximum_coverage $maximum_coverage \
        --max_results $max_results \
        --assembly_quality '${assembly_quality}' \
        --sequences_dir '${sequences_dir}' \
        --outdir $outdir \
        --genome_size_ungapped $genome_size_ungapped

    echo "Sequencing data retrieval complete."
    cat $outdir/sequence_results.json
    """
}
