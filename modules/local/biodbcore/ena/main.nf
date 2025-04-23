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
        path sequence_dir

    output:
        path "ena_results.json", emit: ena_results
        path "**/*.fastq.gz", emit: fastq_files, optional: true
        path "**/*.bam", emit: bam_files, optional: true
        path "**/*.cram", emit: cram_files, optional: true
        path "**/*.sra", emit: sra_files, optional: true
        path "**/*.fast5", emit: fast5_files, optional: true
        path "**/*.pod5", emit: pod5_files, optional: true

    script:
    """
    echo "Checking for existing sequencing data at: $sequence_dir"

    if [ -d "$sequence_dir" ] && find -L "$sequence_dir" -type f \\( \
        -name "*.fastq.gz" -o \
        -name "*.bam" -o \
        -name "*.cram" -o \
        -name "*.sra" \\) | grep -q .; then

        echo "Using local sequencing data from: $sequence_dir for taxonomy ID: $taxonomy_id"

        biodbcore --mode ena \
            --taxonomy_id $taxonomy_id \
            --library_strategy $library_strategy \
            --instrument_platform $instrument_platform \
            --minimum_coverage $minimum_coverage \
            --maximum_coverage $maximum_coverage \
            --max_results $max_results \
            \$( [ -n "$assembly_quality" ] && echo "--assembly_quality $assembly_quality" ) \
            --sequence_dir $sequence_dir \
            --outdir . \
            --genome_size_ungapped $genome_size_ungapped

    else
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
    fi

    echo    "Sequencing data retrieval complete."
    cat     "ena_results.json"
    """
}
