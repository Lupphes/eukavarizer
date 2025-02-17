process REFSEQ_RETRIEVER {
    // Docker
    // container 'eukalyzer/seq_retriver:latest'
    // containerOptions "-v ${params.outdir}/data:/data"

    // Conda
    conda "envs/environment.yml"

    // Settings
    publishDir "$outdir", mode: 'copy', overwrite: false

    input:
        val(taxonomy_id)
        val(outdir)
        val(local_refseq_path)

    output:
        path "$taxonomy_id/refseq_results.json", emit: refseq_json
        path "$taxonomy_id/refseq/*.fna.gz", emit: genome_file
        path "assembly_summary_refseq.parquet", emit: refseq_parquet
        path "assembly_summary_refseq.txt", emit: refseq_summary

    script:
    """
    if [ -f $local_refseq_path ]; then
        echo "Reference genome already exists at: $local_refseq_path for taxonomy ID: $taxonomy_id"

        seq_getter.py --mode refseq --taxonomy_id $taxonomy_id --outdir . --refseq_path $local_refseq_path

    else
        echo "Downloading reference genome for taxonomy ID: $taxonomy_id"

        seq_getter.py --mode refseq --taxonomy_id $taxonomy_id --outdir .
    fi

    cat     './$taxonomy_id/refseq_results.json'
    """
}
