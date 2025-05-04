/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SEQUENCE_PROCESSOR WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    This workflow processes sequencing data and prepares it for structural variant calling:
    1. **BIODBCORE_ENA** – Downloads sequencing data from ENA.
    2. **SAMTOOLS_COLLATEFASTQ** – Converts BAM/CRAM to paired and unpaired FASTQ files.
    3. **QUALITY_CONTROL** – Filters and quality checks the FASTQ files.
    4. **BWAMEM2_MEM** or **BWA_MEM** – Aligns reads to a reference genome.
    6. **SAMTOOLS_INDEX** – Indexes the sorted BAM files.

    Outputs:
    - `fastq_filtered`     – Filtered FASTQ files.
    - `fastq_bam`          – Sorted BAM files (.bam).
    - `fastq_bam_indexes`  – BAM index files (.bai).
    - `multiqc_report`     – MultiQC report file from quality control.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SEQUENCE_FASTQ_CONVERTOR      } from '../../../subworkflows/local/sequence_fastq_convertor/main'
include { QUALITY_CONTROL               } from '../../../subworkflows/local/quality_control/main'
include { SEQUENCE_ALIGNER              } from '../../../subworkflows/local/sequence_aligner/main'
include { SEQUENCE_MERGER               } from '../../../subworkflows/local/sequence_merger/main'

workflow SEQUENCE_PROCESSOR {

    take:
        samplesheet
        reference_genome_unzipped
        reference_genome_bwa_index

    main:
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

    emit:
        fastq_filtered              = QUALITY_CONTROL.out.fastq_filtered    // Filtered FASTQ files
        bam_bai                     = SEQUENCE_MERGER.out.bam_bai           // Sorted, Deduplicated (BAM, BAI) tuples
        multiqc_report              = QUALITY_CONTROL.out.multiqc_report    // Path to MultiQC report
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
