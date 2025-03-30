/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GRIDSS_GRIDSS WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    This workflow calls structural variants (SVs) using GRIDSS for germline
    and somatic samples:
    1. **GRIDSS_GRIDSS** – Detects SVs using local assembly.
    2. **SAMPLE_REHEADER** – Reheaders and renames the output VCF.
    3. **SVYNC** – Synchronizes and refines SV calls.
    4. **GUNZIP** – Decompresses the final VCF file.

    Outputs:
    - `vcf` – Reheaded VCF file
    - `svync_vcf` – Synchronized VCF file
    - `svync_vcfgz` – Gzipped synchronized VCF file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GRIDSS_GRIDSS     } from '../../../modules/nf-core/gridss/gridss/main'
include { SAMPLE_REHEADER   } from '../../../modules/local/sample_regen/main.nf'
include { SVYNC             } from '../../../modules/nf-core/svync/main'
include { GUNZIP            } from '../../../modules/nf-core/gunzip/main'

workflow SV_CALLING_GRIDSS {
    take:
        bam_inputs
        reference_genome_unzipped
        reference_genome_faidx
        reference_genome_bwa_index
        reference_genome_minimap_index

    main:
        name_gridss = "gridss"

        GRIDSS_GRIDSS(
            bam_inputs
                // Don't use on long reads as GRIDSS does not support MINIMAP2 or PacBio/Nanopore raw reads
                .filter { meta, _bam, _bai ->
                    meta.median_bp < params.long_read_threshold
                }
                .map { meta, bam, _bai -> tuple(meta + [id: "${meta.id}_${name_gridss}"], bam) },
            reference_genome_unzipped,
            reference_genome_faidx,
            // reference_genome_minimap_index -> GRIDSS does not support MINIMAP2
            reference_genome_bwa_index
        )

        gridss_result = GRIDSS_GRIDSS.out.vcf

        SAMPLE_REHEADER(
            gridss_result,
            gridss_result.map { meta, _vcf -> "${meta.id}_${name_gridss}" }
        )

        SVYNC(
            SAMPLE_REHEADER.out.vcf
                .join(SAMPLE_REHEADER.out.tbi, by: 0)
                .map { meta, vcf, tbi ->
                    tuple([id: "${meta.id}_svync"], vcf, tbi)
                }
                .combine(
                    Channel.value(file("${projectDir}/assets/svync/${name_gridss}.yaml"))
                )
        )

        GUNZIP(
            SVYNC.out.vcf
        )

    emit:
        vcf = SAMPLE_REHEADER.out.vcf
        vcfgz = SAMPLE_REHEADER.out.vcfgz
        tbi = SAMPLE_REHEADER.out.tbi
        csi = SAMPLE_REHEADER.out.csi
        svync_vcf = GUNZIP.out.gunzip
        svync_vcfgz = SVYNC.out.vcf
        svync_tbi = SVYNC.out.tbi
}
