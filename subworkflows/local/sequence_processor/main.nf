/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SEQUENCE_PROCESSOR WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    This workflow processes sequencing data and prepares it for structural variant calling:
    1. **BIODBCORE_ENA** – Downloads sequencing data from ENA.
    2. **SAMTOOLS_COLLATEFASTQ** – Converts BAM/CRAM to paired and unpaired FASTQ files.
    3. **QUALITY_CONTROL** – Filters and quality checks the FASTQ files.
    4. **BWAMEM2_MEM** or **BWA_MEM** – Aligns reads to a reference genome.
    5. **SAMTOOLS_SORT** – Sorts the BAM files.
    6. **SAMTOOLS_INDEX** – Indexes the sorted BAM files.

    Outputs:
    - `fastq_filtered`     – Filtered FASTQ files.
    - `fastq_bam`          – Sorted BAM files (.bam).
    - `fastq_bam_indexes`  – BAM index files (.bai).
    - `multiqc_report`     – MultiQC report file from quality control.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BIODBCORE_ENA         } from '../../../modules/local/biodbcore/ena/main'

include { GZIP                  } from '../../../modules/local/gzip/main'

include { SRATOOLS_FASTERQDUMP } from '../../../modules/nf-core/sratools/fasterqdump/main'
include { DORADO as DORADO_FAST5    } from '../../../modules/local/dorado/main'
include { DORADO as DORADO_POD5     } from '../../../modules/local/dorado/main'
include { SAMTOOLS_COLLATEFASTQ as BAM_SAMTOOLS_COLLATEFASTQ    } from '../../../modules/nf-core/samtools/collatefastq/main'
include { SAMTOOLS_COLLATEFASTQ as CRAM_SAMTOOLS_COLLATEFASTQ   } from '../../../modules/nf-core/samtools/collatefastq/main'

include { QUALITY_CONTROL       } from '../../../subworkflows/local/quality_control/main'

include { SEQKIT_SIZE           } from '../../../modules/local/seqkit/size/main'
include { MINIMAP2_ALIGN        } from '../../../modules/nf-core/minimap2/align/main'
include { BWAMEM2_MEM           } from '../../../modules/nf-core/bwamem2/mem/main'
include { BWA_MEM               } from '../../../modules/nf-core/bwa/mem/main'

include { SAMTOOLS_MARKDUP      } from '../../../modules/nf-core/samtools/markdup/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_NAMES   } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_COORD   } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX        } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FIXMATE      } from '../../../modules/nf-core/samtools/fixmate/main'

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

        if (sequences_abs_dir != []) {
            // Separate gzipped and uncompressed files
            fastqs_gz = Channel.fromPath("${params.sequence_dir}/**/*.fastq.gz")
            fastqs_unzipped = Channel.fromPath("${params.sequence_dir}/**/*.fastq")

            // Zip the uncompressed ones
            fastqs_zipped = GZIP(
                fastqs_unzipped.map { file -> tuple([id: file.simpleName], file) }
            ).gzip.map { _meta, file -> file }


            // Merge zipped and pre-zipped fastq files
            raw_fastqs = fastqs_gz.mix(fastqs_zipped)
        } else {
            raw_fastqs = BIODBCORE_ENA.out.fastq_files
        }

        raw_sra = (sequences_abs_dir != [] ?
            Channel
                .fromPath("${params.sequence_dir}/**/*.sra")
            : BIODBCORE_ENA.out.sra_files).collect().flatten()

        raw_fast5 = (sequences_abs_dir != [] ?
            Channel
                .fromPath("${params.sequence_dir}/**/*.fast5")
            : BIODBCORE_ENA.out.sra_files).collect().flatten()

        raw_pod5 = (sequences_abs_dir != [] ?
            Channel
                .fromPath("${params.sequence_dir}/**/*.pod5")
            : BIODBCORE_ENA.out.sra_files).collect().flatten()

        raw_bam = (sequences_abs_dir != [] ?
            Channel
                .fromPath("${params.sequence_dir}/**/*.bam")
            : BIODBCORE_ENA.out.bam_files).collect().flatten()

        raw_cram = (sequences_abs_dir != [] ?
            Channel
                .fromPath("${params.sequence_dir}/**/*.cram")
            : BIODBCORE_ENA.out.cram_files).collect().flatten()

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
            .map { file -> tuple([id: file.getParent().getName()], file) }

        grouped_cram = raw_cram
            .flatMap { file -> file instanceof List ? file : [file] }
            .map { file -> tuple([id: file.getParent().getName()], file) }

        grouped_sra = raw_sra
            .flatMap { file -> file instanceof List ? file : [file] }
            .map { file -> tuple([id: file.getParent().getName()], file) }

        grouped_fast5 = raw_fast5
            .flatMap { file -> file instanceof List ? file : [file] }
            .map { file -> tuple([id: file.getParent().getName()], file) }

        grouped_pod5 = raw_pod5
            .flatMap { file -> file instanceof List ? file : [file] }
            .map { file -> tuple([id: file.getParent().getName()], file) }

        SRATOOLS_FASTERQDUMP(
            grouped_sra,
            [],
            []
        )

        DORADO_FAST5(
            grouped_fast5
        )

        DORADO_POD5(
            grouped_pod5
        )

        BAM_SAMTOOLS_COLLATEFASTQ(
            grouped_bam,
            reference_genome_bgzipped.collect(),
            false
        )

        CRAM_SAMTOOLS_COLLATEFASTQ(
            grouped_cram,
            reference_genome_bgzipped.collect(),
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

        // Process SRA parsed reads (paired + unpaired)
        // Paired reads: [[id: sample_id, single_end: true/false], [file1, file2]]
        // Single-end reads: [[id: sample_id, single_end: true/false], [file1]]
        fqdump_reads = SRATOOLS_FASTERQDUMP.out.reads
        .map { meta, fastqs ->
            def single_end = fastqs.size() == 1
            [[id: "${meta.id}", single_end: single_end], fastqs.sort()]
        }

        // Process Dorado parsed reads (fast5) (paired + unpaired)
        // Paired reads: [[id: sample_id, single_end: true/false], [file1, file2]]
        // Single-end reads: [[id: sample_id, single_end: true/false], [file1]]
        dorado_fast5 = DORADO_FAST5.out.fastq
        .map { meta, fastqs ->
            def single_end = fastqs.size() == 1
            [[id: "${meta.id}", single_end: single_end], fastqs ? fastqs : []]
        }

        // Process Dorado parsed reads (pod5) (paired + unpaired)
        // Paired reads: [[id: sample_id, single_end: true/false], [file1, file2]]
        // Single-end reads: [[id: sample_id, single_end: true/false], [file1]]
        dorado_pod5 = DORADO_POD5.out.fastq
        .map { meta, fastqs ->
            def single_end = fastqs.size() == 1
            [[id: "${meta.id}", single_end: single_end], fastqs ? fastqs : []]
        }

        collected_fastqs = grouped_fastqs
            .mix(bam_paired)
            .mix(bam_unpaired)
            .mix(cram_paired)
            .mix(cram_unpaired)
            .mix(fqdump_reads)
            .mix(dorado_fast5)
            .mix(dorado_pod5)

        SEQKIT_SIZE(
            collected_fastqs
        )

        // Add median_bp to the metadata
        tagged_collected_fastqs = SEQKIT_SIZE.out.median_bp.map { meta, fastq, median_bp ->
            def length = median_bp.text.trim().toInteger()
            tuple(meta + [median_bp: length], fastq)
        }

        tagged_collected_fastqs.view()

        QUALITY_CONTROL(
            tagged_collected_fastqs
        )

        minimap2_bam = QUALITY_CONTROL.out.fastq_filtered
            .filter { meta, _fastq ->
                params.minimap2_flag && (meta.median_bp > params.long_read_threshold)
            }

        bwa_bam = QUALITY_CONTROL.out.fastq_filtered
            .filter { meta, _fastq ->
                !params.minimap2_flag || (meta.median_bp <= params.long_read_threshold)
            }

        MINIMAP2_ALIGN(
            minimap2_bam,
            reference_genome_unzipped.collect(),
            true,
            [],
            false,
            false
        )

        sorted_indexed_bed = params.bwamem2 ?
            BWAMEM2_MEM(bwa_bam, reference_genome_bwa_index.collect(), reference_genome_unzipped.collect(), false) :
            BWA_MEM(bwa_bam, reference_genome_bwa_index.collect(), reference_genome_unzipped.collect(), false)

        mixed_bam_inputs = sorted_indexed_bed.bam
            .mix(MINIMAP2_ALIGN.out.bam)
            .map { meta, bam ->
                tuple(meta + [id: "${meta.id}_mixed"], bam)
            }

        // Sort paired BAM files
        SAMTOOLS_SORT_NAMES(
            mixed_bam_inputs.filter { meta, _fastq -> !meta.single_end },
            reference_genome_unzipped.collect()
        )

        // FIXMATE for paired BAM files
        SAMTOOLS_FIXMATE(
            SAMTOOLS_SORT_NAMES.out.bam,
        )

        // Sort mixed BAM files
        SAMTOOLS_SORT_COORD(
            SAMTOOLS_FIXMATE.out.bam.map { meta, bam ->
                tuple(meta + [id: meta.id.replaceFirst(/_mixed$/, '_coo')], bam)
            }.mix(
                mixed_bam_inputs.filter { meta, _fastq -> meta.single_end }
            ),
            reference_genome_unzipped.collect()
        )

        // Mark duplicates of paired BAM files
        mark_dup = SAMTOOLS_MARKDUP(
            SAMTOOLS_SORT_COORD.out.bam.map { meta, bam ->
                tuple(meta + [id: meta.id.replaceFirst(/_coo$/, '_mar')], bam)
            }.filter { meta, _fastq -> !meta.single_end },
            reference_genome_unzipped.collect()
        )

        // Mix sorted BAM files with unpaired reads with marked duplicates
        samtools_result = mark_dup.bam.map { meta, bam ->
            tuple(meta + [id: meta.id.replaceFirst(/_mar$/, '_sas')], bam)
        }
        .mix(SAMTOOLS_SORT_COORD.out.bam.filter { meta, _fastq -> meta.single_end })

        // Sort and index BAM files
        SAMTOOLS_INDEX(
            samtools_result
        )

    emit:
        fastq_filtered              = QUALITY_CONTROL.out.fastq_filtered    // Filtered FASTQ files
        fastq_bam                   = samtools_result                       // Sorted BAM files (.bam)
        fastq_bam_indexes           = SAMTOOLS_INDEX.out.bai                // BAM index files (.bai)
        multiqc_report              = QUALITY_CONTROL.out.multiqc_report    // Path to MultiQC report
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
