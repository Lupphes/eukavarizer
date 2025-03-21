/*
 * GRIDSS_GRIDSS
 *   Advanced, local assembly-based SV caller for short-read data
 *   (both germline and somatic).
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

    main:
        name_gridss = "gridss"

        GRIDSS_GRIDSS(
            bam_inputs.map { meta, bam, _bai -> tuple(meta + [id: "${meta.id}_${name_gridss}"], bam) },
            reference_genome_unzipped,
            reference_genome_faidx,
            reference_genome_bwa_index
        )

        SAMPLE_REHEADER(
            GRIDSS_GRIDSS.out.vcf,
            GRIDSS_GRIDSS.out.vcf.map { meta, _vcf -> "${meta.id}_${name_gridss}" }
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
