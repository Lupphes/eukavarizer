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
    val(remove_headers)

    output:
    tuple val(meta), path("reheadered_${new_name}.vcf"), emit: vcf, optional: true
    tuple val(meta), path("reheadered_${new_name}.vcf.gz"), emit: vcfgz, optional: true
    tuple val(meta), path("reheadered_${new_name}.vcf.gz.tbi"), emit: tbi, optional: true
    tuple val(meta), path("reheadered_${new_name}.vcf.gz.csi"), emit: csi, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # Determine file extension and ensure compatibility
    if [[ "${vcf}" == *.vcf.gz ]]; then
        VCF_INPUT="${vcf}"
    else
        VCF_INPUT="input.vcf.gz"
        bgzip -c "${vcf}" > "\${VCF_INPUT}"
        tabix -p vcf "\${VCF_INPUT}"
    fi

    # Check if the VCF is readable (must have at least a #CHROM line)
    if bcftools view "\${VCF_INPUT}" | grep -q '^#CHROM'; then
        # Get sample names count
        SAMPLE_COUNT=\$(bcftools query -l "\${VCF_INPUT}" | wc -l)

        if [[ "\${SAMPLE_COUNT}" -gt 0 ]]; then
            if [[ "${remove_headers}" == "true" ]]; then
                # Strip/filter FILTER headers with empty or malformed descriptions and sanitize
                bcftools view -h "\${VCF_INPUT}" | awk '
                    /^##FILTER=/ {
                        match(\$0, /Description="[^"]+"/)
                        if (RSTART > 0) {
                            desc = substr(\$0, RSTART + 12, RLENGTH - 13)
                            # Mask greater-than symbol
                            gsub(/>/, "\\\\&gt;", desc)
                            sub(/Description="[^"]+"/, "Description=\\"" desc "\\"")
                        }
                    }
                    { print }
                ' > filtered_header.hdr
            else
                bcftools view -h "\${VCF_INPUT}" > filtered_header.hdr
            fi

            # Extract original sample names and prepend new_name
            bcftools query -l "\${VCF_INPUT}" | awk -v prefix="[${new_name}]" '{print prefix \$0}' > sample_names_${new_name}.txt

            # Create a new header with updated sample names
            bcftools reheader -h filtered_header.hdr -s sample_names_${new_name}.txt -o reheadered_${new_name}.vcf.gz "\${VCF_INPUT}"

            # Index the new VCF (both TBI and CSI)
            tabix -p vcf reheadered_${new_name}.vcf.gz
            bcftools index --csi reheadered_${new_name}.vcf.gz

            # Also provide an uncompressed VCF
            bcftools view reheadered_${new_name}.vcf.gz -Ov -o reheadered_${new_name}.vcf
        else
            echo "No samples found in VCF. Skipping reheadering."
        fi
    else
        echo "Invalid or empty VCF. Skipping reheadering."
    fi

    # Capture versions
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -n 1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """


    stub:
    """
    touch reheadered_${new_name}.vcf.gz
    touch reheadered_${new_name}.vcf.gz.tbi
    touch reheadered_${new_name}.vcf.gz.csi
    touch versions.yml
    """
}
