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
    - `vcf`           – Reheaded VCF file.
    - `vcfgz`         – Gzipped reheaded VCF file.
    - `tbi`           – Tabix index (.tbi) for the reheaded VCF.
    - `csi`           – CSI index (.csi) for the reheaded VCF (if generated).
    - `svync_vcf`     – Decompressed synchronized VCF file.
    - `svync_vcfgz`   – Gzipped synchronized VCF file.
    - `svync_tbi`     – Tabix index for the synchronized VCF.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GRIDSS_GRIDSS     } from '../../../modules/nf-core/gridss/gridss/main'
include { GRIDSS_ANNOTATE   } from '../../../modules/local/gridss/annotate/main.nf'
include { SVYNC             } from '../../../modules/nf-core/svync/main'
include { SAMPLE_REHEADER   } from '../../../modules/local/sample_regen/main.nf'

workflow SV_CALLING_GRIDSS {
    take:
        bam_inputs
        reference_genome_unzipped
        reference_genome_faidx
        reference_genome_bwa_index

    main:
        name_gridss = "gridss"

        GRIDSS_GRIDSS(
            bam_inputs
                // Don't use on long reads as GRIDSS does not support MINIMAP2 or PacBio/Nanopore raw reads
                .filter { meta, _bam, _bai ->
                    meta.median_bp < params.long_read_threshold
                }
                .map { meta, bam, _bai -> tuple(meta + [id: "${meta.id}-${name_gridss}"], bam) },
            reference_genome_unzipped,
            reference_genome_faidx,
            // reference_genome_minimap_index -> GRIDSS does not support MINIMAP2
            reference_genome_bwa_index
        )

        if (params.gridss_annotate) {
            gridss_annotate = GRIDSS_ANNOTATE(
                GRIDSS_GRIDSS.out.vcf
            ).vcf
        } else {
            gridss_annotate = GRIDSS_GRIDSS.out.vcf
        }

        SVYNC(
            gridss_annotate
                .map { meta, vcf ->
                    tuple(meta + [id: "${meta.id}-svync"], vcf, [])
                }
                .combine(
                    Channel.value(file("${projectDir}/assets/svync/${name_gridss}.yaml"))
                )
        )

        SAMPLE_REHEADER(
            SVYNC.out.vcf,
            name_gridss,
            ""
        )

    emit:
        vcf = SAMPLE_REHEADER.out.vcf
        vcfgz =  SAMPLE_REHEADER.out.vcfgz
        tbi = SAMPLE_REHEADER.out.tbi
        csi = SAMPLE_REHEADER.out.csi
}
