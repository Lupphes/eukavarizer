process SAMPLE_REHEADER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0' :
        'biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf)
    val(algorithm)
    // TODO: Make this a parameter
    val(suffix)

    output:
    tuple val(meta), path("re-${meta.sample}-${meta.platform}-${meta.single_end}-${algorithm}${suffix}.vcf"), emit: vcf, optional: true
    tuple val(meta), path("re-${meta.sample}-${meta.platform}-${meta.single_end}-${algorithm}${suffix}.vcf.gz"), emit: vcfgz, optional: true
    tuple val(meta), path("re-${meta.sample}-${meta.platform}-${meta.single_end}-${algorithm}${suffix}.vcf.gz.tbi"), emit: tbi, optional: true
    tuple val(meta), path("re-${meta.sample}-${meta.platform}-${meta.single_end}-${algorithm}${suffix}.vcf.gz.csi"), emit: csi, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

script:
def algorihm_sample_id = "${meta.sample}-${meta.platform}-${meta.single_end}-${algorithm}${suffix}"
def prefix_algorihm_sample_id = "re-${algorihm_sample_id}"
"""
# Ensure input is compressed and indexed
if [[ "${vcf}" == *.vcf.gz ]]; then
    VCF_INPUT="${vcf}"
else
    VCF_INPUT="input.vcf.gz"
    bgzip --threads ${task.cpus} -c "${vcf}" > "\${VCF_INPUT}"
    tabix --threads ${task.cpus} -p vcf "\${VCF_INPUT}"
fi

# Count samples
SAMPLE_COUNT=\$(bcftools query -l "\${VCF_INPUT}" | wc -l)

if [[ "\${SAMPLE_COUNT}" -eq 1 ]]; then
    # Rename sample
    echo "${meta.sample}" > new_sample.txt
    bcftools reheader -s new_sample.txt -o "renamed-\${VCF_INPUT}" "\${VCF_INPUT}"

    # Inject INFO header lines
    bcftools view -h "renamed-\${VCF_INPUT}" | grep -E "^##" > tmp.vcf
    echo '##INFO=<ID=EUK_CALLER,Number=1,Type=String,Description="SV calling algorithm used">' >> tmp.vcf
    echo '##INFO=<ID=EUK_PLATFORM,Number=1,Type=String,Description="Sequencing platform">' >> tmp.vcf
    echo '##INFO=<ID=EUK_SE,Number=1,Type=String,Description="Single-end or paired-end">' >> tmp.vcf

    # Add column header
    bcftools view -h "renamed-\${VCF_INPUT}" | grep -E "^#CHROM" >> tmp.vcf

    # Add INFO values to each variant
    bcftools view -H "renamed-\${VCF_INPUT}" | \\
    awk -v algo="${algorithm}" -v plat="${meta.platform}" -v se="${meta.single_end}" \\
    'BEGIN {FS=OFS="\t"} {
        if (\$8 == ".") {
            \$8 = "EUK_CALLER=" algo ";EUK_PLATFORM=" plat ";EUK_SE=" se
        } else {
            \$8 = \$8 ";EUK_CALLER=" algo ";EUK_PLATFORM=" plat ";EUK_SE=" se
        }
        print
    }' >> tmp.vcf

    # Compress and index output
    bgzip --threads ${task.cpus} -c tmp.vcf > ${prefix_algorihm_sample_id}.vcf.gz
    tabix --threads ${task.cpus} ${prefix_algorihm_sample_id}.vcf.gz
    bcftools index --csi ${prefix_algorihm_sample_id}.vcf.gz

    # Uncompressed version
    bcftools view ${prefix_algorihm_sample_id}.vcf.gz -Ov -o ${prefix_algorihm_sample_id}.vcf

elif [[ "\${SAMPLE_COUNT}" -eq 0 ]]; then
    echo "No samples found in VCF. Skipping reheadering."

else
    echo "Error: VCF contains more than one sample (\${SAMPLE_COUNT}). Cannot proceed." >&2
    exit 1
fi

# Capture versions
cat <<EOF > versions.yml
"${task.process}":
  bcftools: \$(bcftools --version | head -n 1 | sed 's/^.*bcftools //')
EOF
"""
    stub:
    def algorihm_sample_id = "${meta.sample}"
    def prefix_algorihm_sample_id = "reheadered-${algorihm_sample_id}"
    """
    touch ${prefix_algorihm_sample_id}.vcf.gz
    touch ${prefix_algorihm_sample_id}.vcf.gz.tbi
    touch ${prefix_algorihm_sample_id}.vcf.gz.csi
    touch versions.yml
    """
}
