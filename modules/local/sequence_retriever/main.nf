process SEQUENCE_RETRIEVER {
    // Docker
    // container 'eukalyzer/seq_retriver:latest'
    // containerOptions "-v ${params.outdir}/data:/data"

    // Conda
    conda "envs/environment.yml"

    // Settings
    publishDir "$outdir", mode: 'copy'

    input:
        val(taxonomy_id)
        val(genome_size_ungapped)
        val(outdir)
        val(library_strategy)
        val(instrument_platform)
        val(minimum_coverage)
        val(maximum_coverage)
        val(max_results)
        val(assembly_quality)
        val(sequence_dir)

    output:
        path "$taxonomy_id/ena_results.json", emit: ena_results
        path "$taxonomy_id/sequences/*/*.fastq.gz", emit: sequence_files


    script:
    """
    echo "Checking for existing sequencing data at: $sequence_dir"

    if [ -d "$sequence_dir" ] && find "$sequence_dir" -type f -name "*.fastq.gz" | grep -q .; then
        echo "Using local sequencing data from: $sequence_dir for taxonomy ID: $taxonomy_id"

        mkdir -p "$taxonomy_id/sequences"
        cp -r "${sequence_dir}/" "$taxonomy_id/"

        seq_getter.py --mode ena \
            --taxonomy_id $taxonomy_id \
            --library_strategy $library_strategy \
            --instrument_platform $instrument_platform \
            --minimum_coverage $minimum_coverage \
            --maximum_coverage $maximum_coverage \
            --max_results $max_results \
            \$( [ -n "$assembly_quality" ] && echo "--assembly_quality $assembly_quality" ) \
            --sequences_dir $sequence_dir \
            --outdir . \
            --genome_size_ungapped $genome_size_ungapped

    else
        echo "No local sequencing data found. Downloading sequencing reads for taxonomy ID: $taxonomy_id"

        seq_getter.py --mode ena \
            --taxonomy_id $taxonomy_id \
            --library_strategy $library_strategy \
            --instrument_platform $instrument_platform \
            --minimum_coverage $minimum_coverage \
            --maximum_coverage $maximum_coverage \
            --max_results $max_results \
            \$( [ -n "$assembly_quality" ] && echo "--assembly_quality $assembly_quality" ) \
            --outdir . \
            --genome_size_ungapped $genome_size_ungapped
    fi

    echo    "Sequencing data retrieval complete."
    cat     $taxonomy_id/ena_results.json
    """
}
