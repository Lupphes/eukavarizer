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
    - `vcf` – Reheaded VCF file
    - `svync_vcf` – Synchronized VCF file
    - `svync_vcfgz` – Gzipped synchronized VCF file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CUTESV            } from '../../../modules/nf-core/cutesv/main'
include { SAMPLE_REHEADER   } from '../../../modules/local/sample_regen/main.nf'
include { SVYNC             } from '../../../modules/nf-core/svync/main'
include { GUNZIP            } from '../../../modules/nf-core/gunzip/main'

workflow SV_CALLING_CUTESV {
    take:
        bam_inputs
        reference_genome_bgzipped

    main:
        name_cutesv = "cutesv"

        CUTESV(
            bam_inputs.map { meta, bam, bai ->
                tuple(meta + [id: "${meta.id}_${name_cutesv}"], bam, bai)
            },
            reference_genome_bgzipped.map { meta, fasta ->
                tuple(meta + [id: "${meta.id}_${name_cutesv}"], fasta)
            }
        )

        SAMPLE_REHEADER(
            CUTESV.out.vcf,
            CUTESV.out.vcf.map { meta, _vcf -> "${meta.id}_${name_cutesv}" }
        )

        SVYNC(
            SAMPLE_REHEADER.out.vcf
                .join(SAMPLE_REHEADER.out.tbi, by: 0)
                .map { meta, vcf, tbi ->
                    tuple([id: "${meta.id}_svync"], vcf, tbi)
                }
                .combine(
                    Channel.value(file("${projectDir}/assets/svync/${name_cutesv}.yaml"))
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
