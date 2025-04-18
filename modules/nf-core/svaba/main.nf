
process SVABA {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/svaba:1.1.0--h7d7f7ad_2':
        'biocontainers/svaba:1.1.0--h7d7f7ad_2' }"

    input:
    tuple val(meta), path(tumorbam), path(tumorbai), path(normalbam), path(normalbai)
    tuple val(meta1), path(fasta)
    tuple val(meta2), path(fasta_fai)
    tuple val(meta3), path(bwa_index)
    tuple val(meta4), path(dbsnp)
    tuple val(meta5), path(dbsnp_tbi)
    tuple val(meta6), path(regions)

    output:
    tuple val(meta), path("*.svaba.sv.vcf.gz")                        , emit: sv, optional: true
    tuple val(meta), path("*.svaba.indel.vcf.gz")                     , emit: indel, optional: true
    tuple val(meta), path("*.svaba.germline.indel.vcf.gz")            , emit: germ_indel, optional: true
    tuple val(meta), path("*.svaba.germline.sv.vcf.gz")               , emit: germ_sv, optional: true
    tuple val(meta), path("*.svaba.somatic.indel.vcf.gz")             , emit: som_indel,  optional: true
    tuple val(meta), path("*.svaba.somatic.sv.vcf.gz")                , emit: som_sv, optional: true
    tuple val(meta), path("*.svaba.unfiltered.sv.vcf.gz")             , emit: unfiltered_sv, optional: true
    tuple val(meta), path("*.svaba.unfiltered.indel.vcf.gz")          , emit: unfiltered_indel, optional: true
    tuple val(meta), path("*.svaba.unfiltered.germline.indel.vcf.gz") , emit: unfiltered_germ_indel, optional: true
    tuple val(meta), path("*.svaba.unfiltered.germline.sv.vcf.gz")    , emit: unfiltered_germ_sv, optional: true
    tuple val(meta), path("*.svaba.unfiltered.somatic.indel.vcf.gz")  , emit: unfiltered_som_indel,  optional: true
    tuple val(meta), path("*.svaba.unfiltered.somatic.sv.vcf.gz")     , emit: unfiltered_som_sv, optional: true
    tuple val(meta), path("*.bps.txt.gz")                             , emit: raw_calls
    tuple val(meta), path("*.discordants.txt.gz")                     , emit: discordants, optional: true
    tuple val(meta), path("*.log")                                    , emit: log
    path "versions.yml"                                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def bamlist = normalbam ? "-t ${tumorbam} -n ${normalbam}" : "-t ${tumorbam}"
    dbsnp   = dbsnp ? "--dbsnp-vcf ${dbsnp}" : ""
    regions = regions ? "--region ${regions}" : ""

    """
    ref_file=\$(basename "${fasta}")
    ref_prefix=\${ref_file%.*}

    if [ ! -f "./\${ref_file}" ]; then
        cp "${fasta}" .
    fi

    for ext in amb ann bwt pac sa; do
        cp "${bwa_index}/\${ref_prefix}.\$ext" "\${ref_file}.\$ext"
    done

    svaba \\
        run \\
        $bamlist \\
        --threads $task.cpus \\
        $dbsnp \\
        --id-string ${prefix} \\
        --reference-genome \${ref_file} \\
        --g-zip \\
        $regions \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svaba: \$(echo \$(svaba --version 2>&1) | sed 's/[^0-9.]*\\([0-9.]*\\).*/\\1/' )
    END_VERSIONS
    """
    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bps.txt.gz
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svaba: \$(echo \$(svaba --version 2>&1) | sed 's/[^0-9.]*\\([0-9.]*\\).*/\\1/' )
    END_VERSIONS
    """
}
