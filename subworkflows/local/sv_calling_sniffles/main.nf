/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SNIFFLES WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    This workflow detects structural variants (SVs) from long-read data using Sniffles:
    1. **SNIFFLES** – Calls SVs using split-read and coverage-based methods.
    2. **SAMPLE_REHEADER** – Renames and reheaders the VCF output.
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

include { SNIFFLES          } from '../../../modules/nf-core/sniffles/main'
include { SVYNC             } from '../../../modules/nf-core/svync/main'
include { SAMPLE_REHEADER   } from '../../../modules/local/sample_regen/main.nf'

workflow SV_CALLING_SNIFFLES {
    take:
        bam_inputs
        reference_genome_bgzipped

    main:
        name_sniffles = "sniffles"

        SNIFFLES(
            bam_inputs
            .map { meta, bam, bai ->
                tuple(meta + [id: "${meta.id}-${name_sniffles}"], bam, bai)
            },
            reference_genome_bgzipped,
            [[],[]],
            true,
            false
        )

        SVYNC(
            SNIFFLES.out.vcf
                .join(SNIFFLES.out.tbi, by: 0)
                .map { meta, vcf, tbi ->
                    tuple(meta + [id: "${meta.id}-${name_sniffles}-svync"], vcf, tbi)
                }
                .combine(
                    Channel.value(file("${projectDir}/assets/svync/${name_sniffles}.yaml"))
                )
        )

        SAMPLE_REHEADER(
            SVYNC.out.vcf,
            name_sniffles,
            ""
        )

    emit:
        vcf = SAMPLE_REHEADER.out.vcf
        vcfgz =  SAMPLE_REHEADER.out.vcfgz
        tbi = SAMPLE_REHEADER.out.tbi
        csi = SAMPLE_REHEADER.out.csi
}
