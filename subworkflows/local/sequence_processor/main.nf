/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BIODBCORE_REFSEQ      } from '../../../modules/local/biodbcore/refseq/main'
include { BIODBCORE_ENA         } from '../../../modules/local/biodbcore/ena/main'

include { SAMTOOLS_COLLATEFASTQ as BAM_SAMTOOLS_COLLATEFASTQ    } from '../../../modules/nf-core/samtools/collatefastq/main'
include { SAMTOOLS_COLLATEFASTQ as CRAM_SAMTOOLS_COLLATEFASTQ   } from '../../../modules/nf-core/samtools/collatefastq/main'

include { QUALITY_CONTROL       } from '../../../subworkflows/local/quality_control/main'

include { BWA_MEM               } from '../../../modules/nf-core/bwa/mem/main'

include { SAMTOOLS_SORT         } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX        } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FAIDX        } from '../../../modules/nf-core/samtools/faidx/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SEQUENCE_PROCESSOR {

    take:
        taxonomy_id
        outdir
        sequences_abs_dir
        reference_genome_ungapped_size
        reference_genome_unzipped
        reference_genome_bgzipped
        reference_genome_bwa_index

    main:
        ch_library_strategy    = params.library_strategy ? Channel.value(params.library_strategy) : Channel.value([])
        ch_instrument_platform = params.instrument_platform ? Channel.value(params.instrument_platform) : Channel.value([])
        ch_minimum_coverage    = params.minimum_coverage ? Channel.value(params.minimum_coverage) : Channel.value(1)
        ch_maximum_coverage    = params.maximum_coverage ? Channel.value(params.maximum_coverage) : Channel.value("")
        ch_max_results         = params.max_results ? Channel.value(params.max_results) : Channel.value(1)
        ch_assembly_quality    = params.assembly_quality ? Channel.value(params.assembly_quality) : Channel.value("")

        BIODBCORE_ENA(
            taxonomy_id,
            reference_genome_ungapped_size,
            outdir,
            ch_library_strategy.map     { it.join(' ').trim() },
            ch_instrument_platform.map  { it.join(' ').trim() },
            ch_minimum_coverage,
            ch_maximum_coverage,
            ch_max_results,
            ch_assembly_quality,
            sequences_abs_dir
        )

        raw_fastqs  = BIODBCORE_ENA.out.fastq_files
        raw_bam     = BIODBCORE_ENA.out.bam_files
        raw_cram    = BIODBCORE_ENA.out.cram_files


        // Paired reads: [[id: sample_id, single_end: true/false], [file1, file2]]
        // Single-end reads: [[id: sample_id, single_end: true/false], [file1]]
        grouped_fastqs = raw_fastqs
            .flatMap { file ->
                file instanceof List ? file : [file]
            }
            .map { file -> tuple(file.getParent().getName(), file) }
            .groupTuple()
            .map { run_accession, files ->
                def paired = files.findAll { f ->
                    def base = f.toString().replaceAll(/_[12]\.fastq\.gz$/, '')
                    files.any { it.toString() == "${base}_1.fastq.gz" } &&
                    files.any { it.toString() == "${base}_2.fastq.gz" }
                }

                def is_paired = paired.size() == 2

                if (is_paired) {
                    // Paired-end
                    tuple([id: run_accession, single_end: !is_paired], paired.sort())
                } else {
                    // Single-end
                    def unpaired = files - paired
                    tuple([id: run_accession, single_end: !is_paired], unpaired.sort())
                }
            }

        grouped_bam = raw_bam
            .flatMap { file -> file instanceof List ? file : [file] }
            .map { file -> tuple([id: file.getParent().getName()], file.toString()) }
            .groupTuple()

        grouped_cram = raw_cram
            .flatMap { file -> file instanceof List ? file : [file] }
            .map { file -> tuple([id: file.getParent().getName()], file.toString()) }
            .groupTuple()


        BAM_SAMTOOLS_COLLATEFASTQ(
            grouped_bam,
            reference_genome_bgzipped,
            false
        )

        CRAM_SAMTOOLS_COLLATEFASTQ(
            grouped_cram,
            reference_genome_bgzipped,
            false
        )

        // Process BAM reads (paired + unpaired)
        // Paired reads: [[id: sample_id, single_end: true/false], [file1, file2]]
        // Single-end reads: [[id: sample_id, single_end: true/false], [file1]]
        bam_paired = BAM_SAMTOOLS_COLLATEFASTQ.out.fastq
            .filter { _meta, files -> files.size() == 2 }
            .map { meta, fastqs ->
                [[id: "${meta.id}_paired", single_end: false], fastqs.sort()]
            }

        bam_unpaired = BAM_SAMTOOLS_COLLATEFASTQ.out.fastq_singleton
            .map { meta, fastqs ->
                [[id: "${meta.id}_unpaired", single_end: true], fastqs ? fastqs : []]
            }

        // Process CRAM reads (paired + unpaired)
        // Paired reads: [[id: sample_id, single_end: true/false], [file1, file2]]
        // Single-end reads: [[id: sample_id, single_end: true/false], [file1]]
        cram_paired = CRAM_SAMTOOLS_COLLATEFASTQ.out.fastq
            .filter { _meta, files -> files.size() == 2 }
            .map { meta, fastqs ->
                [[id: "${meta.id}_paired", single_end: false], fastqs.sort()]
            }

        cram_unpaired = CRAM_SAMTOOLS_COLLATEFASTQ.out.fastq_singleton
            .map { meta, fastqs ->
                [[id: "${meta.id}_unpaired", single_end: true], fastqs ? fastqs : []]
            }

        collected_fastqs = grouped_fastqs
            .mix(bam_paired)
            .mix(bam_unpaired)
            .mix(cram_paired)
            .mix(cram_unpaired)

        QUALITY_CONTROL(
            collected_fastqs
        )

        BWA_MEM(
            QUALITY_CONTROL.out.fastq_filtered,
            reference_genome_bwa_index,
            reference_genome_unzipped,
            true
        )

        SAMTOOLS_SORT(
            BWA_MEM.out.bam.map { meta, bam ->
                tuple(meta + [id: "${meta.id}_sas"], bam)
            },
            reference_genome_unzipped
        )

        SAMTOOLS_INDEX(
            SAMTOOLS_SORT.out.bam
        )


    emit:
        fastq_filtered              = QUALITY_CONTROL.out.fastq_filtered
        fastq_bam                   = SAMTOOLS_SORT.out.bam             // Sorted BAM files (.bam)
        fastq_bam_indexes           = SAMTOOLS_INDEX.out.bai            // BAM index files (.bai)
        // fasta_file       = genome_unzipped.map { it[1] }    // Unzipped FASTA file
        // fasta_index      = ch_fasta_index                   // FASTA index (.fai) for unzipped FASTA
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
