/*
 * DYSGU
 *   Machine-learning approach for structural variant calling from short reads.
 */
include { DYSGU             } from '../../../modules/nf-core/dysgu/main'
include { SAMPLE_REHEADER   } from '../../../modules/local/sample_regen/main.nf'
include { SVYNC             } from '../../../modules/nf-core/svync/main'
include { GUNZIP            } from '../../../modules/nf-core/gunzip/main'

workflow SV_CALLING_DYSGU {
    take:
        bam_inputs
        reference_genome_unzipped
        reference_genome_faidx

    main:
        name_dysgu = "dysgu"

        first = bam_inputs
            .map { meta, bam, bai ->
                tuple(meta + [id: "${meta.id}_${name_dysgu}"], bam, bai)
            }

        second = reference_genome_unzipped
            .join(reference_genome_faidx, by: 0)
            .map { meta, fasta, fai ->
                tuple(meta + [id: "${meta.id}_${name_dysgu}"], fasta, fai)
            }

        DYSGU(
            first,
            second,
        )

        SAMPLE_REHEADER(
            DYSGU.out.vcf,
            DYSGU.out.vcf.map { meta, _vcf -> "${meta.id}_${name_dysgu}" }
        )

        SVYNC(
            SAMPLE_REHEADER.out.vcf
                .join(SAMPLE_REHEADER.out.tbi, by: 0)
                .map { meta, vcf, tbi ->
                    tuple([id: "${meta.id}_svync"], vcf, tbi)
                }
                .combine(
                    Channel.value(file("${projectDir}/assets/svync/${name_dysgu}.yaml"))
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
