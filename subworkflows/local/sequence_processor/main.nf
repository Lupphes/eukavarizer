/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: SEQUENCE_PROCESSOR
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Complete sequencing data processing pipeline from raw data to analysis-ready BAM files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Description:
        This subworkflow orchestrates the complete processing of sequencing data, from
        format conversion through quality control, alignment, and post-processing. It
        produces analysis-ready BAM files with quality metrics for downstream variant calling.

    Processing Steps:
        1. Convert input formats to FASTQ (SEQUENCE_FASTQ_CONVERTOR)
        2. Perform quality control and filtering (QUALITY_CONTROL)
        3. Align reads to reference genome (SEQUENCE_ALIGNER)
        4. Merge, deduplicate, and recalibrate BAM files (SEQUENCE_MERGER)
        5. Generate alignment statistics (SAMTOOLS_STATS)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Input Channels (take):
        samplesheet               tuple(meta, files)     Input sequencing data
        reference_genome_unzipped tuple(meta, fasta)     Uncompressed reference genome
        reference_genome_bwa_index tuple(meta, index)    BWA/BWA-MEM2 index files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Output Channels (emit):
        fastq_filtered            tuple(meta, fastq)     Filtered FASTQ files
        bam_bai                   tuple(meta, bam, bai)  Sorted, deduplicated BAM with index
        samtools_stats            tuple(meta, stats)     Samtools alignment statistics
        multiqc_report            path(html)             MultiQC quality control report
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Author:   Ondřej Sloup (Lupphes)
    Contact:  ondrej.sloup@protonmail.com
    GitHub:   @Lupphes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SEQUENCE_FASTQ_CONVERTOR      } from '../../../subworkflows/local/sequence_fastq_convertor/main'
include { QUALITY_CONTROL               } from '../../../subworkflows/local/quality_control/main'
include { SEQUENCE_ALIGNER              } from '../../../subworkflows/local/sequence_aligner/main'
include { SEQUENCE_MERGER               } from '../../../subworkflows/local/sequence_merger/main'
include { SAMTOOLS_STATS                } from '../../../modules/nf-core/samtools/stats/main'

workflow SEQUENCE_PROCESSOR {

    take:
        samplesheet
        reference_genome_unzipped
        reference_genome_bwa_index

    main:
        ch_versions = Channel.empty()

        SEQUENCE_FASTQ_CONVERTOR(
            samplesheet,
            reference_genome_unzipped,
            reference_genome_bwa_index
        )

        QUALITY_CONTROL(
            SEQUENCE_FASTQ_CONVERTOR.out.tagged_collected_fastqs
        )

        SEQUENCE_ALIGNER(
            QUALITY_CONTROL.out.fastq_filtered,
            reference_genome_unzipped,
            reference_genome_bwa_index
        )

        SEQUENCE_MERGER(
            SEQUENCE_ALIGNER.out.bam_mapped,
            reference_genome_unzipped,
            reference_genome_bwa_index
        )

        SAMTOOLS_STATS(
            SEQUENCE_MERGER.out.bam_bai,
            reference_genome_unzipped
        )

        ch_versions = ch_versions.mix(SEQUENCE_FASTQ_CONVERTOR.out.versions)
        ch_versions = ch_versions.mix(QUALITY_CONTROL.out.versions)
        ch_versions = ch_versions.mix(SEQUENCE_ALIGNER.out.versions)
        ch_versions = ch_versions.mix(SEQUENCE_MERGER.out.versions)
        ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions.first())

    emit:
        fastq_filtered              = QUALITY_CONTROL.out.fastq_filtered    // Filtered FASTQ files
        bam_bai                     = SEQUENCE_MERGER.out.bam_bai           // Sorted, Deduplicated (BAM, BAI) tuples
        samtools_stats              = SAMTOOLS_STATS.out.stats              // Samtools stats files
        multiqc_report              = QUALITY_CONTROL.out.multiqc_report    // Path to MultiQC report
        versions                    = ch_versions                           // Software versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
