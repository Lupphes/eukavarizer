/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TIDDIT_SV WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    This workflow detects structural variants (SVs) using TIDDIT and processes the output:
    1. **TIDDIT_SV** – Detects SVs from BAM files using a coverage-based approach.
    2. **SAMPLE_REHEADER** – Renames and reheaders the VCF output.
    3. **SVYNC** – Synchronizes and refines SV calls.
    4. **GUNZIP** – Decompresses the final VCF file.

    Outputs:
    - `vcf` – Reheaded VCF file
    - `svync_vcf` – Synchronized VCF file
    - `svync_vcfgz` – Gzipped synchronized VCF file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { TIDDIT_SV         } from '../../../modules/nf-core/tiddit/sv/main'
include { SAMPLE_REHEADER   } from '../../../modules/local/sample_regen/main.nf'
include { SVYNC             } from '../../../modules/nf-core/svync/main'
include { GUNZIP            } from '../../../modules/nf-core/gunzip/main'


workflow SV_CALLING_TIDDIT {
    take:
        bam_inputs
        reference_genome_bgzipped
        reference_genome_bwa_index

    main:
        name_tiddit = "tiddit"

        TIDDIT_SV(
            bam_inputs.map { meta, bam, bai ->
                tuple(meta + [id: "${meta.id}_${name_tiddit}"], bam, bai)
            },
            reference_genome_bgzipped,
            reference_genome_bwa_index.map { meta, bwa ->
                tuple(meta + [id: "${meta.id}_${name_tiddit}"], bwa)
            }
        )

        SAMPLE_REHEADER(
            TIDDIT_SV.out.vcf,
            TIDDIT_SV.out.vcf.map { meta, _vcf -> "${meta.id}_${name_tiddit}" }
        )

        SVYNC(
            SAMPLE_REHEADER.out.vcf
                .join(SAMPLE_REHEADER.out.tbi, by: 0)
                .map { meta, vcf, tbi ->
                    tuple([id: "${meta.id}_svync"], vcf, tbi)
                }
                .combine(
                    Channel.value(file("${projectDir}/assets/svync/${name_tiddit}.yaml"))
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
