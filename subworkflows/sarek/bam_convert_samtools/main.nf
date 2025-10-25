//
// BAM/CRAM to FASTQ conversion, supports both single-end and paired-end data
// This workflow is adapted from the Sarek pipeline with modifications to support single-end data.
// Original Sarek version: https://raw.githubusercontent.com/nf-core/sarek/refs/heads/master/subworkflows/local/bam_convert_samtools/main.nf
//
// MODIFICATIONS FROM UPSTREAM:
// - Added conditional logic to handle single-end data (line 60-72)
// - Original Sarek version only supported paired-end data
// - Used updated module versions
//

include { SAMTOOLS_VIEW         as SAMTOOLS_VIEW_MAP_MAP     } from '../../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_VIEW         as SAMTOOLS_VIEW_UNMAP_UNMAP } from '../../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_VIEW         as SAMTOOLS_VIEW_UNMAP_MAP   } from '../../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_VIEW         as SAMTOOLS_VIEW_MAP_UNMAP   } from '../../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_MERGE        as SAMTOOLS_MERGE_UNMAP      } from '../../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_INDEX                                     } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_COLLATEFASTQ as COLLATE_FASTQ_UNMAP       } from '../../../modules/nf-core/samtools/collatefastq/main'
include { SAMTOOLS_COLLATEFASTQ as COLLATE_FASTQ_MAP         } from '../../../modules/nf-core/samtools/collatefastq/main'
include { CAT_FASTQ                                          } from '../../../modules/nf-core/cat/fastq/main'

workflow BAM_CONVERT_SAMTOOLS {
    take:
    input       // channel: [meta, alignment (BAM or CRAM), index (optional)]
    fasta       // optional: reference file if CRAM format and reference not in header
    fasta_fai
    interleaved // value: true/false

    main:
    versions = channel.empty()

    SAMTOOLS_INDEX(input)
    bam_bai = input.join(SAMTOOLS_INDEX.out.bai.mix(SAMTOOLS_INDEX.out.crai), by: 0)

    // MAP - MAP
    SAMTOOLS_VIEW_MAP_MAP(bam_bai, fasta, [], [])

    // UNMAP - UNMAP
    SAMTOOLS_VIEW_UNMAP_UNMAP(bam_bai, fasta, [], [])

    // UNMAP - MAP
    SAMTOOLS_VIEW_UNMAP_MAP(bam_bai, fasta, [], [])

    // MAP - UNMAP
    SAMTOOLS_VIEW_MAP_UNMAP(bam_bai, fasta, [], [])

    // Merge UNMAP
    all_unmapped_bam = SAMTOOLS_VIEW_UNMAP_UNMAP.out.bam.mix(SAMTOOLS_VIEW_UNMAP_UNMAP.out.cram)
        .join(SAMTOOLS_VIEW_UNMAP_MAP.out.bam.mix(SAMTOOLS_VIEW_UNMAP_MAP.out.cram), failOnDuplicate: true, remainder: true)
        .join(SAMTOOLS_VIEW_MAP_UNMAP.out.bam.mix(SAMTOOLS_VIEW_MAP_UNMAP.out.cram), failOnDuplicate: true, remainder: true)
        .map{ meta, unmap_unmap, unmap_map, map_unmap -> [ meta, [ unmap_unmap, unmap_map, map_unmap ] ] }

    SAMTOOLS_MERGE_UNMAP(all_unmapped_bam, fasta, fasta_fai, [[],[]])

    // Collate & convert unmapped
    COLLATE_FASTQ_UNMAP(SAMTOOLS_MERGE_UNMAP.out.bam.mix(SAMTOOLS_MERGE_UNMAP.out.cram), fasta, interleaved)

    // Collate & convert mapped
    COLLATE_FASTQ_MAP(SAMTOOLS_VIEW_MAP_MAP.out.bam.mix(SAMTOOLS_VIEW_MAP_MAP.out.cram), fasta, interleaved)

    // Join Mapped & unmapped fastq
    reads_to_concat = COLLATE_FASTQ_MAP.out.fastq
        .join(COLLATE_FASTQ_UNMAP.out.fastq, failOnDuplicate: true, failOnMismatch: true)
        .map{ meta, mapped_reads, unmapped_reads ->
            // Handle both single-end and paired-end data
            // COLLATE_FASTQ outputs: paired-end = [R1, R2], single-end = single file
            // CAT_FASTQ expects: paired-end = [file1_R1, file1_R2, file2_R1, file2_R2, ...]
            //                   single-end = [file1, file2, file3, ...]

            if (meta.single_end) {
                // Single-end: Flatten arrays into single list
                def mapped = mapped_reads instanceof List ? mapped_reads : [mapped_reads]
                def unmapped = unmapped_reads instanceof List ? unmapped_reads : [unmapped_reads]
                [ meta, (mapped + unmapped).flatten() ]
            } else {
                // Paired-end: Interleave R1 and R2 files properly
                // Ensure we have arrays (defensive programming)
                def mapped = mapped_reads instanceof List ? mapped_reads : [mapped_reads]
                def unmapped = unmapped_reads instanceof List ? unmapped_reads : [unmapped_reads]

                // For paired-end, we expect [R1, R2] from each COLLATE_FASTQ
                // Output should be: [mapped_R1, mapped_R2, unmapped_R1, unmapped_R2]
                [ meta, [ mapped[0], mapped[1], unmapped[0], unmapped[1] ] ]
            }
        }

    // Concatenate Mapped_R1 with Unmapped_R1 and Mapped_R2 with Unmapped_R2
    CAT_FASTQ(reads_to_concat)
    reads = CAT_FASTQ.out.reads

    // Gather versions of all tools used
    versions = versions.mix(CAT_FASTQ.out.versions)
    versions = versions.mix(COLLATE_FASTQ_MAP.out.versions)
    versions = versions.mix(COLLATE_FASTQ_UNMAP.out.versions)
    versions = versions.mix(SAMTOOLS_MERGE_UNMAP.out.versions)
    versions = versions.mix(SAMTOOLS_VIEW_MAP_MAP.out.versions)
    versions = versions.mix(SAMTOOLS_VIEW_MAP_UNMAP.out.versions)
    versions = versions.mix(SAMTOOLS_VIEW_UNMAP_MAP.out.versions)
    versions = versions.mix(SAMTOOLS_VIEW_UNMAP_UNMAP.out.versions)

    emit:
    reads

    versions
}
