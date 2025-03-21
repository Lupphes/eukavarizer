/*
 * SVABA
 *   Local assembly-based SV and indel detection (single-sample or tumor/normal).
 */
include { SVABA             } from '../../../modules/nf-core/svaba/main'
include { SAMPLE_REHEADER   } from '../../../modules/local/sample_regen/main.nf'
include { SVYNC             } from '../../../modules/nf-core/svync/main'
include { GUNZIP            } from '../../../modules/nf-core/gunzip/main'

workflow SV_CALLING_SVABA {
    take:
        bam_inputs
        reference_genome_unzipped
        reference_genome_faidx
        reference_genome_bwa_index

    main:
        name_svaba = "svaba"

        first = bam_inputs
            .map { meta, bam, bai ->
                tuple(meta + [id: "${meta.id}_svaba"], bam, bai, [], [])
            }

        SVABA(
            first,
            reference_genome_unzipped,
            reference_genome_faidx,
            reference_genome_bwa_index,
            [[],[]],
            [[],[]],
            [[],[]],
        )

        SAMPLE_REHEADER(
            SVABA.out.germ_sv,
            SVABA.out.germ_sv.map { meta, _vcf -> "${meta.id}_${name_svaba}" }
        )

        SVYNC(
            SAMPLE_REHEADER.out.vcf
                .join(SAMPLE_REHEADER.out.tbi, by: 0)
                .map { meta, vcf, tbi ->
                    tuple([id: "${meta.id}_svync"], vcf, tbi)
                }
                .combine(
                    Channel.value(file("${projectDir}/assets/svync/${name_svaba}.yaml"))
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
