process DOWNLOADER {
    container 'eukalyzer/downloader:latest'
    errorStrategy 'terminate'
    debug true

    publishDir "${params.outdir}/sequences", mode: 'copy'

    // Mounts output directory inside the container
    containerOptions "-v ${params.outdir}/data:/data"

    input:
    val(user_input)

    output:
    path "sequences/*.fastq.gz", emit: fastq_files

    script:
    """
    echo "Downloading sequencing data with the following parameters:"
    echo "Taxonomy ID: ${params.taxonomy_id}"
    echo "Library Strategies: ${params.library_strategy.join(' ')}"
    echo "Instrument Platforms: ${params.instrument_platform.join(' ')}"
    echo "Coverage Range: ${params.minimum_coverage} - ${params.maximum_coverage}"
    echo "Max Results: ${params.max_results}"
    echo "Reference Genome: ${params.ref_genome_path}"
    echo "Sequences Directory: ${params.sequences_dir}"

    python -u /app/main.py \
        --taxonomy_id ${params.taxonomy_id} \
        --library_strategy '${params.library_strategy.join(",")}' \
        --instrument_platform '${params.instrument_platform.join(",")}' \
        --minimum_coverage ${params.minimum_coverage} \
        --maximum_coverage ${params.maximum_coverage} \
        --max_results ${params.max_results} \
        --ref_genome_path ${params.ref_genome_path} \
        --sequences_dir ${params.sequences_dir}

    echo "Download complete."

    mkdir -p sequences
    cp -v ${params.sequences_dir}/*.fastq.gz sequences/

    echo "Files copied to Nextflow work directory:"
    ls -lah sequences/
    """
}
