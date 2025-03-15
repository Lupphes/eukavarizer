process REFSEQ_RETRIEVER {
    tag "$taxonomy_id"
    conda "${moduleDir}/environment.yml"

    // Settings
    publishDir "$outdir", mode: 'copy', overwrite: false

    input:
        val(taxonomy_id)
        val(outdir)
        val(genome_file)

    output:
        path "$taxonomy_id/refseq_results.json", emit: refseq_json
        path "$taxonomy_id/refseq/*.fna.gz", emit: genome_file
        path "assembly_summary_refseq.parquet", emit: refseq_parquet, optional: true
        path "assembly_summary_refseq.txt", emit: refseq_summary, optional: true

    script:
    """
    if [ -s "${genome_file}" ]; then
        echo "Reference genome already exists at: $genome_file for taxonomy ID: $taxonomy_id"

        mkdir -p "$taxonomy_id/refseq"
        cp "${genome_file}" "$taxonomy_id/refseq/"

        seq_getter.py --mode refseq --taxonomy_id $taxonomy_id --outdir . --genome_file $genome_file

    else
        echo "Downloading reference genome for taxonomy ID: $taxonomy_id"

        seq_getter.py --mode refseq --taxonomy_id $taxonomy_id --outdir .
    fi

    cat     './$taxonomy_id/refseq_results.json'
    """
}
