process SAMPLE_REHEADER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0' :
        'biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf)
    val(new_name)

    output:
    tuple val(meta), path("reheadered_${new_name}.vcf"), emit: vcf
    tuple val(meta), path("reheadered_${new_name}.vcf.gz"), emit: vcfgz
    tuple val(meta), path("reheadered_${new_name}.vcf.gz.tbi"), emit: tbi
    tuple val(meta), path("reheadered_${new_name}.vcf.gz.csi"), emit: csi

    script:
    """
    # Determine file extension and ensure compatibility
    if [[ "${vcf}" == *.vcf.gz ]]; then
        VCF_INPUT="${vcf}"
    else
        VCF_INPUT="input.vcf.gz"
        bgzip -c "${vcf}" > "\${VCF_INPUT}"
        tabix -p vcf "\${VCF_INPUT}"
    fi

    # Get sample names count
    SAMPLE_COUNT=\$(bcftools query -l "\${VCF_INPUT}" | wc -l)

    if [[ "\${SAMPLE_COUNT}" -gt 0 ]]; then
        # Extract original sample names and prepend new_name
        bcftools query -l "\${VCF_INPUT}" | awk -v prefix="[${new_name}]" '{print prefix \$0}' > sample_names_${new_name}.txt

        # Create a new header with updated sample names
        bcftools reheader -s sample_names_${new_name}.txt -o reheadered_${new_name}.vcf.gz "\${VCF_INPUT}"

        # Index the new VCF (both TBI and CSI)
        tabix -p vcf reheadered_${new_name}.vcf.gz  # Generates .tbi
        bcftools index --csi reheadered_${new_name}.vcf.gz  # Generates .csi

        # Also provide an uncompressed VCF
        bcftools view reheadered_${new_name}.vcf.gz -Ov -o reheadered_${new_name}.vcf
    else
        echo "No samples found in VCF. Skipping reheadering."
        cp "\${VCF_INPUT}" "reheadered_${new_name}.vcf.gz"
        tabix -p vcf "reheadered_${new_name}.vcf.gz"
        bcftools index --csi "reheadered_${new_name}.vcf.gz"
        bcftools view "reheadered_${new_name}.vcf.gz" -Ov -o "reheadered_${new_name}.vcf"
    fi
    """
}
