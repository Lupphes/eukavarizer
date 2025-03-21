/*
 * SNIFFLES
 *   Long-read (PacBio/ONT) structural variant caller using split-read logic.
 */
include { SNIFFLES          } from '../../../modules/nf-core/sniffles/main'
include { SAMPLE_REHEADER   } from '../../../modules/local/sample_regen/main.nf'
include { SVYNC             } from '../../../modules/nf-core/svync/main'
include { GUNZIP            } from '../../../modules/nf-core/gunzip/main'

workflow SV_CALLING_SNIFFLES {
    take:
        bam_inputs
        reference_genome_bgzipped

    main:
        name_sniffles = "sniffles"

        SNIFFLES(
            bam_inputs
            .map { meta, bam, bai ->
                tuple(meta + [id: "${meta.id}_${name_sniffles}"], bam, bai)
            },
            reference_genome_bgzipped.map { meta, fasta ->
                tuple(meta + [id: "${meta.id}_${name_sniffles}"], fasta)
            },
            [[],[]],
            true,
            false
        )

        SAMPLE_REHEADER(
            SNIFFLES.out.vcf,
            SNIFFLES.out.vcf.map { meta, _bfc -> "${meta.id}" }
        )

        SVYNC(
            SAMPLE_REHEADER.out.vcf
                .join(SAMPLE_REHEADER.out.tbi, by: 0)
                .map { meta, vcf, tbi ->
                    tuple([id: "${meta.id}_svync"], vcf, tbi)
                }
                .combine(
                    Channel.value(file("${projectDir}/assets/svync/${name_sniffles}.yaml"))
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
