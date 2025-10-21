process SAMPLE_REHEADER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/47/474a5ea8dc03366b04df884d89aeacc4f8e6d1ad92266888e7a8e7958d07cde8/data':
        'community.wave.seqera.io/library/bcftools_htslib:0a3fa2654b52006f' }"

    input:
    tuple val(meta), path(vcf)
    val(algorithm)

    output:
    tuple val(meta), path("${prefix}.vcf"), emit: vcf, optional: true
    tuple val(meta), path("${prefix}.vcf.gz"), emit: vcfgz, optional: true
    tuple val(meta), path("${prefix}.vcf.gz.tbi"), emit: tbi, optional: true
    tuple val(meta), path("${prefix}.vcf.gz.csi"), emit: csi, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def suffix = task.ext.suffix ?: ''
    prefix = task.ext.prefix ?: "re-${meta.sample}-${meta.platform}-${meta.single_end}-${algorithm}${suffix}"
    def meta_platform = meta.platform ?: "unknown"
    """
    #!/bin/bash
    set -euo pipefail

    # Ensure input is compressed and indexed
    if [[ "${vcf}" == *.vcf.gz ]]; then
        VCF_INPUT="${vcf}"
    else
        VCF_INPUT="input.vcf.gz"
        bgzip --threads ${task.cpus} -c "${vcf}" > "\${VCF_INPUT}"
        tabix --threads ${task.cpus} -p vcf "\${VCF_INPUT}"
    fi

    # Count samples (use || true to handle empty VCFs in strict mode)
    SAMPLE_COUNT=\$(bcftools query -l "\${VCF_INPUT}" 2>/dev/null | wc -l || true)

    if [[ "\${SAMPLE_COUNT}" -eq 1 ]]; then
        # Rename sample using bcftools reheader
        echo "${meta.sample}" > new_sample.txt
        bcftools reheader -s new_sample.txt "\${VCF_INPUT}" -o temp_reheadered.vcf.gz

        # Stream processing - add INFO headers and values
        {
            # Output all header lines except #CHROM line
            bcftools view --no-version -h temp_reheadered.vcf.gz | grep -v "^#CHROM"

            # Add custom INFO headers
            cat <<'INFO_HEADERS'
##INFO=<ID=EUK_CALLER,Number=1,Type=String,Description="SV calling algorithm used">
##INFO=<ID=EUK_PLATFORM,Number=1,Type=String,Description="Sequencing platform">
##INFO=<ID=EUK_SE,Number=1,Type=String,Description="Single-end or paired-end">
INFO_HEADERS

            # Output #CHROM line
            bcftools view --no-version -h temp_reheadered.vcf.gz | grep "^#CHROM"

            # Add INFO values to each variant
            bcftools view --no-version -H temp_reheadered.vcf.gz | awk -v algo="${algorithm}" -v plat="${meta_platform}" -v se="${meta.single_end}" \\
                'BEGIN {FS=OFS="\\t"} {
                    info_add = "EUK_CALLER=" algo ";EUK_PLATFORM=" plat ";EUK_SE=" se
                    if (\$8 == "." || \$8 == "") {
                        \$8 = info_add
                    } else {
                        \$8 = \$8 ";" info_add
                    }
                    print
                }'
        } | bgzip --threads ${task.cpus} > ${prefix}.vcf.gz

        # Index the output
        tabix --threads ${task.cpus} -p vcf ${prefix}.vcf.gz
        bcftools index --threads ${task.cpus} --csi ${prefix}.vcf.gz

        # Create uncompressed version
        bcftools view --no-version -Ov -o ${prefix}.vcf ${prefix}.vcf.gz

        # Cleanup
        rm -f temp_reheadered.vcf.gz new_sample.txt

    elif [[ "\${SAMPLE_COUNT}" -eq 0 ]]; then
        echo "WARNING: No samples found in VCF (sites-only). Skipping reheadering." >&2

    else
        echo "ERROR: VCF contains \${SAMPLE_COUNT} samples (expected 1). Cannot proceed." >&2
        exit 1
    fi

    # Capture versions
    cat <<EOF > versions.yml
"${task.process}":
    bcftools: \$(bcftools --version 2>&1 | head -n 1 | sed 's/^.*bcftools //')
    tabix: \$(tabix --version 2>&1 | head -n 1 | sed 's/^.*tabix (htslib) //' | sed 's/ .*//')
EOF
    """

    stub:
    def suffix = task.ext.suffix ?: ''
    prefix = task.ext.prefix ?: "re-${meta.sample}-${meta.platform}-${meta.single_end}-${algorithm}${suffix}"
    """
    touch ${prefix}.vcf
    echo -n | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    touch ${prefix}.vcf.gz.csi
    touch versions.yml
    """
}
