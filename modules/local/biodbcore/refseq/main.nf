process BIODBCORE_REFSEQ {
    tag "$taxonomy_id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "luppo/biodbcore:latest"

    input:
        val  taxonomy_id
        val  outdir
        path reference_genome

    output:
        path "refseq_results.json", emit: json
        path "*.{fna.gz,fa.gz}", emit: reference_genome, optional: true
        path "assembly_summary_refseq.parquet", emit: refseq_parquet, optional: true
        path "assembly_summary_refseq.txt", emit: refseq_summary, optional: true

    script:
        if (reference_genome) {
        """
            echo "Reference genome already exists at: $reference_genome for taxonomy ID: $taxonomy_id"

            biodbcore --mode refseq --taxonomy_id $taxonomy_id --outdir . --reference_genome $reference_genome

            cat 'refseq_results.json'
        """
        } else {
        """
            echo "Downloading reference genome for taxonomy ID: $taxonomy_id"

            biodbcore --mode refseq --taxonomy_id $taxonomy_id --outdir .
            cat 'refseq_results.json'
        """
        }
}
