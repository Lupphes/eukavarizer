include { MINIMAP2_ALIGN        } from '../../../modules/nf-core/minimap2/align/main'
include { BWAMEM2_MEM           } from '../../../modules/nf-core/bwamem2/mem/main'
include { BWA_MEM               } from '../../../modules/nf-core/bwa/mem/main'

workflow SEQUENCE_ALIGNER {
    take:
        fastq_filtered
        reference_genome_unzipped
        reference_genome_bwa_index

    main:
        // THIS HAS BEEN ADAPTED FROM NF_CORE/SAREK
        // STEP 1: MAPPING READS TO REFERENCE GENOME
        // First, we must calculate number of lanes for each sample (meta.n_fastq)
        // This is needed to group reads from the same sample together using groupKey to avoid stalling the workflow
        // when reads from different samples are mixed together
        fastq_filtered.map { meta, reads ->
                [ meta.subMap('patient', 'sample', 'sex', 'status', 'platform', 'single_end'), reads ]
            }
            .groupTuple()
            .map { meta, reads ->
                meta + [ n_fastq: reads.size() ] // We can drop the FASTQ files now that we know how many there are
            }
            .set { reads_grouping_key }

        reads_for_alignment = fastq_filtered.map{ meta, reads ->
            // Update meta.id to meta.sample no multiple lanes or splitted fastqs
            if (meta.size * meta.num_lanes == 1) [ meta + [ id:meta.sample ], reads ]
            else [ meta, reads ]
        }

        minimap2_bam = reads_for_alignment
            .filter { meta, _fastq ->
                params.minimap2_flag && (
                    (meta.platform && meta.platform == 'ont' || meta.platform == 'pacbio') ||
                    (!meta.platform && meta.median_bp > params.long_read_threshold)
                )
            }

        bwa_bam = reads_for_alignment
            .filter { meta, _fastq ->
                !params.minimap2_flag || (
                    (meta.platform && meta.platform == 'illumina') ||
                    (!meta.platform && meta.median_bp <= params.long_read_threshold)
                )
            }

        MINIMAP2_ALIGN(
            minimap2_bam,
            reference_genome_unzipped.collect(),
            true,
            "",
            false,
            true
        )

        sorted_indexed_bed = params.bwamem2 ?
            BWAMEM2_MEM(bwa_bam, reference_genome_bwa_index.collect(), reference_genome_unzipped.collect(), true) :
            BWA_MEM(bwa_bam, reference_genome_bwa_index.collect(), reference_genome_unzipped.collect(), true)

        mixed_bam_inputs = sorted_indexed_bed.bam.mix(MINIMAP2_ALIGN.out.bam)

        // Grouping the bams from the same samples not to stall the workflow
        // Use groupKey to make sure that the correct group can advance as soon as it is complete
        // and not stall the workflow until all reads from all channels are mapped
        bam_mapped = mixed_bam_inputs
            .combine(reads_grouping_key) // Creates a tuple of [ meta, bam, reads_grouping_key ]
            .filter {
                meta1, _bam, meta2 -> meta1.sample == meta2.sample
                && meta1.single_end == meta2.single_end
                && meta1.platform == meta2.platform
            }
            // Add n_fastq and other variables to meta
            .map { meta1, bam, meta2 ->
                [ meta1 + meta2, bam ]
            }
            // Manipulate meta map to remove old fields and add new ones
            .map { meta, bam ->
                [ meta - meta.subMap('id', 'read_group', 'data_type', 'num_lanes', 'read_group', 'size', 'median_bp') + [ data_type: 'bam', id: "${meta.sample}-${meta.platform}-${meta.single_end}" ], bam ]
            }
            // Create groupKey from meta map
            .map { meta, bam ->
                [ groupKey( meta, meta.n_fastq), bam ]
            }
            // Group
            .groupTuple()

    emit:
        bam_mapped = bam_mapped
}
