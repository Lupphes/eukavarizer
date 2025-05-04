/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CUTESV WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    This workflow calls structural variants (SVs) using CuteSV for long-read data:
    1. **CUTESV** – Detects SVs using split-read and coverage-based methods.
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

include { CUTESV            } from '../../../modules/nf-core/cutesv/main'
include { SVYNC             } from '../../../modules/nf-core/svync/main'
include { SAMPLE_REHEADER   } from '../../../modules/local/sample_regen/main.nf'

workflow SV_CALLING_CUTESV {
    take:
        bam_inputs
        reference_genome_bgzipped

    main:
        name_cutesv = "cutesv"

        CUTESV(
            bam_inputs.map { meta, bam, bai ->
                tuple(meta + [id: "${meta.id}-${name_cutesv}"], bam, bai)
            },
            reference_genome_bgzipped
        )

        SVYNC(
            CUTESV.out.vcf
                .map { meta, vcf ->
                    tuple(meta + [id: "${meta.id}-svync"], vcf, [])
                }
                .combine(
                    Channel.value(file("${projectDir}/assets/svync/${name_cutesv}.yaml"))
                )
        )

        SAMPLE_REHEADER(
            SVYNC.out.vcf,
            name_cutesv,
            ""
        )

    emit:
        vcf = SAMPLE_REHEADER.out.vcf
        vcfgz =  SAMPLE_REHEADER.out.vcfgz
        tbi = SAMPLE_REHEADER.out.tbi
        csi = SAMPLE_REHEADER.out.csi
}
