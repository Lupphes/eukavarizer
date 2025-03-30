/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DELLY_CALL WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    This workflow calls structural variants (SVs) using Delly:
    1. **DELLY_CALL** – Detects deletions, duplications, inversions, and translocations.
    2. **SAMPLE_REHEADER** – Reheaders and renames the output VCF.
    3. **SVYNC** – Synchronizes and refines SV calls.
    4. **GUNZIP** – Decompresses the final VCF file.

    Outputs:
    - `vcf` – Reheaded VCF file
    - `svync_vcf` – Synchronized VCF file
    - `svync_vcfgz` – Gzipped synchronized VCF file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DELLY_CALL        } from '../../../modules/nf-core/delly/call/main.nf'
include { SAMPLE_REHEADER   } from '../../../modules/local/sample_regen/main.nf'
include { SVYNC             } from '../../../modules/nf-core/svync/main'
include { GUNZIP            } from '../../../modules/nf-core/gunzip/main'

workflow SV_CALLING_DELLY {

    take:
        fastq_bam
        fastq_bam_indexes
        reference_genome_bgzipped
        reference_genome_faidx

    main:
        name_delly = "delly"

        ch_delly_input = fastq_bam
            .join(fastq_bam_indexes, by: 0)
            .map { meta, bam, bai ->
                tuple(meta + [id: "${meta.id}_${name_delly}"], bam, bai, [], [], [])
            }

        DELLY_CALL(
            ch_delly_input,
            reference_genome_bgzipped,
            reference_genome_faidx
        )

        SAMPLE_REHEADER(
            DELLY_CALL.out.bcf, // vcf
            DELLY_CALL.out.bcf.map { meta, _bcf -> "${meta.id}" }
        )

        SVYNC(
            SAMPLE_REHEADER.out.vcf
                .join(SAMPLE_REHEADER.out.tbi, by: 0)
                .map { meta, vcf, tbi ->
                    tuple([id: "${meta.id}_svync"], vcf, tbi)
                }
                .combine(
                    Channel.value(file("${projectDir}/assets/svync/${name_delly}.yaml"))
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
