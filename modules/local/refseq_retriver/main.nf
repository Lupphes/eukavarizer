    process REFSEQ_RETRIEVER {
        // Docker
        // container 'eukalyzer/seq_retriver:latest'
        // containerOptions "-v ${params.outdir}/data:/data"

        // Conda
        conda "envs/environment.yml"

        // Settings
        errorStrategy 'terminate'
        publishDir "${params.outdir}/refseq", mode: 'copy'

        input:
            val(taxonomy_id)
            val(outdir)
            val(local_refseq_path)

    output:
        path "$outdir/refseq_results.json", emit: refseq_json

    script:
    """
    mkdir -p refseq

    if [ -f $local_refseq_path ]; then
        echo "Reference genome already exists at: $local_refseq_path for taxonomy ID: $taxonomy_id"

        python /app/main.py --mode refseq --taxonomy_id $taxonomy_id --outdir $outdir --refseq_path $local_refseq_path

        cp -v $local_refseq_path refseq/

    else
        echo "Downloading reference genome for taxonomy ID: $taxonomy_id"

        python /app/main.py --mode refseq --taxonomy_id $taxonomy_id --outdir $outdir

        cp -v $outdir/*.fna.gz refseq/
    fi

    cat $outdir/refseq_results.json
    """
}
