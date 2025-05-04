include { SAMTOOLS_INDEX as INDEX_MERGE_BAM } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as INDEX_MARKDUPLICATES } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_MERGE                    } from '../../../modules/nf-core/samtools/merge/main'

include { GATK4SPARK_MARKDUPLICATES         } from '../../../modules/nf-core/gatk4spark/markduplicates/main'
include { GATK4SPARK_BASERECALIBRATOR       } from '../../../modules/nf-core/gatk4spark/baserecalibrator/main'
include { GATK4SPARK_APPLYBQSR              } from '../../../modules/nf-core/gatk4spark/applybqsr/main'

include { GATK4_MARKDUPLICATES              } from '../../../modules/nf-core/gatk4/markduplicates/main'
include { GATK4_BASERECALIBRATOR            } from '../../../modules/nf-core/gatk4/baserecalibrator/main'
include { GATK4_APPLYBQSR                   } from '../../../modules/nf-core/gatk4/applybqsr/main'

workflow SEQUENCE_MERGER {
    take:
        bam_mapped
        reference_genome_unzipped
        reference_genome_bwa_index

    main:
        if (params.deduplicate_flag) {
            // TODO: Explore why SPARK did not work
            GATK4_MARKDUPLICATES(
                bam_mapped,
                reference_genome_unzipped.map{ _meta, fasta -> [ fasta ] }.collect(),
                reference_genome_bwa_index.map{ _meta, fasta_fai -> [ fasta_fai ] }.collect()
            )

            INDEX_MARKDUPLICATES(GATK4_MARKDUPLICATES.out.bam)

            bam_bai = GATK4_MARKDUPLICATES.out.bam.join(INDEX_MARKDUPLICATES.out.bai, failOnDuplicate: true, failOnMismatch: true)
        } else {
            bam_to_merge = bam_mapped.branch{ meta, bam ->
                // bam is a list, so use bam.size() to asses number of intervals
                single:   bam.size() <= 1
                    return [ meta, bam[0] ]
                multiple: bam.size() > 1
            }

            // Only when using intervals
            SAMTOOLS_MERGE(bam_to_merge.multiple, [ [ id:'null' ], []], [ [ id:'null' ], []])

            // Mix intervals and no_intervals channels together
            bam_all = SAMTOOLS_MERGE.out.bam.mix(bam_to_merge.single)

            // Index bam
            INDEX_MERGE_BAM(bam_all)

            bam_bai = bam_all.join(INDEX_MERGE_BAM.out.bai, failOnDuplicate: true, failOnMismatch: true)
        }

        if (params.recalibrate_flag && params.known_sites && params.known_sites_tbi) {
            GATK4_BASERECALIBRATOR(
                bam_bai.map{ meta, bam, bai -> tuple(meta, bam, bai, []) },
                reference_genome_unzipped.collect(),
                reference_genome_bwa_index.collect(),
                [ [ id:'dict' ], [] ],
                params.known_sites,
                params.known_sites_tbi
            )

            bam_bai_table = bam_bai.join(GATK4_BASERECALIBRATOR.out.table, failOnDuplicate: true, failOnMismatch: true)

            GATK4_APPLYBQSR(
                bam_bai_table.map{ meta, bam, bai, table -> tuple(meta, bam, bai, table, []) },
                reference_genome_unzipped.map{ _meta, fasta -> [ fasta ] }.collect(),
                reference_genome_bwa_index.map{ _meta, fasta_fai -> [ fasta_fai ] }.collect(),
                [ [ id:'dict' ], [] ],
            )

            re_bam_bai = GATK4_APPLYBQSR.out.bam
        } else {
            re_bam_bai = bam_bai
        }
    emit:
        bam_bai = re_bam_bai

}
