/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES FOR ALIGNMENT & INDEXING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GUNZIP         } from '../modules/nf-core/gunzip/main'
include { BWA_INDEX      } from '../modules/nf-core/bwa/index/main'
include { BWA_MEM        } from '../modules/nf-core/bwa/mem/main'
include { SAMTOOLS_SORT  } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FAIDX } from '../modules/nf-core/samtools/faidx/main'
include { TABIX_BGZIP    } from '../modules/nf-core/tabix/bgzip/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW: INPUT_GENERATION
    - Aligns reads to a reference genome
    - Sorts and indexes BAM files
    - Prepares reference genome for downstream analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow INPUT_GENERATION {

    take:
        grouped_fastqs  // FASTQ input files (Paired and Unpaired)
        ch_genome_file  // Reference genome (FASTA.GZ)

    main:
        log.info "🚀 Starting INPUT_GENERATION workflow"

        // Debug: Display input files
        grouped_fastqs.view { "DEBUG: Received grouped FASTQ files -> ${it}" }
        ch_genome_file.view { "DEBUG: Received genome file -> ${it}" }

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            STEP 1: Decompress Reference Genome (if needed)
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        ch_unzipped_fasta = GUNZIP(
            ch_genome_file.map { file -> tuple([id: file.simpleName.replaceFirst(/\.gz$/, '')], file) }
        ).gunzip

        ch_unzipped_fasta.view { "DEBUG: Decompressed genome -> ${it}" }

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            STEP 2: Generate BWA Index
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        ch_bwa_index = BWA_INDEX(ch_unzipped_fasta).index

        ch_bwa_index.view { "DEBUG: BWA index -> ${it}" }

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            STEP 3: Prepare FASTQ Files for Alignment
            - Groups paired-end reads together
            - Assigns sample IDs
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        // Debug: View the raw grouped fastqs before processing
        grouped_fastqs.view { "DEBUG: Raw grouped FASTQs -> ${it}" }

        ch_prepared_fastqs = grouped_fastqs.map { run_accession, paired, unpaired ->
            def meta = [id: run_accession, single_end: paired.size() != 2] // Ensure correct meta map
            def reads = paired.size() == 2 ? [paired[0], paired[1]] : [paired[0]] // Handle single/paired

            return tuple(meta, reads)
        }

        // Debug: Check structure before sending to BWA_MEM
        ch_prepared_fastqs.view { "DEBUG: Formatted FASTQ tuples -> ${it}" }

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            STEP 4: Align Reads to Reference Genome (BWA-MEM)
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        // Debug: View BWA_MEM inputs before execution
        ch_prepared_fastqs.view { "DEBUG: BWA_MEM input -> ${it}" }
        ch_bwa_index.view { "DEBUG: BWA index input -> ${it}" }
        ch_unzipped_fasta.view { "DEBUG: Unzipped FASTA input -> ${it}" }

        ch_aligned_bam = BWA_MEM(
            ch_prepared_fastqs,  // This now contains meta + list of reads
            ch_bwa_index.map { meta, index -> tuple(meta, index) },
            ch_unzipped_fasta.map { meta, fasta -> tuple(meta, fasta) },
            true
        ).bam

        // Debug: View output from BWA_MEM
        ch_aligned_bam.view { "DEBUG: BWA_MEM output BAM -> ${it}" }


        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            STEP 5: Sort BAM Files
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        ch_sorted_bam = SAMTOOLS_SORT(
            ch_aligned_bam.map { meta, bam ->
                def new_meta = meta.clone()
                new_meta.id = "${meta.id}_sorted"
                new_meta.prefix = "${meta.id}_sorted"
                tuple(new_meta, bam)
            },
            ch_unzipped_fasta.map { meta, fasta -> tuple(meta, fasta) }
        ).bam

        ch_sorted_bam.view { "DEBUG: Sorted BAM files -> ${it}" }

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            STEP 6: Index BAM Files
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        ch_bam_index = SAMTOOLS_INDEX(ch_sorted_bam).bai

        ch_bam_index.view { "DEBUG: BAM index files -> ${it}" }

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            STEP 7: Generate FASTA Index (.fai)
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        ch_fasta_index = SAMTOOLS_FAIDX(
            ch_unzipped_fasta.map { tuple(it[0], it[1]) }, // (meta, fasta)
            [[], []]
        ).fai

        ch_fasta_index.view { "DEBUG: FASTA index (.fai) -> ${it}" }

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            STEP 8: Compress FASTA with bgzip
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        ch_bgzipped_fasta_tuple = TABIX_BGZIP(
            ch_unzipped_fasta.map { meta, fasta -> tuple(meta, fasta) }
        )

        ch_bgzipped_fasta = ch_bgzipped_fasta_tuple.output.map { meta, file -> file }
        ch_gzi_index = ch_bgzipped_fasta_tuple.gzi

        ch_bgzipped_fasta.view { "DEBUG: BGZIPPED FASTA -> ${it}" }
        ch_gzi_index.view { "DEBUG: BGZIP INDEX FILE -> ${it}" }

    emit:
        bam_files        = ch_sorted_bam
        bam_indexes      = ch_bam_index
        fasta_file       = ch_unzipped_fasta.map { meta, fasta -> fasta }
        fasta_index      = ch_fasta_index
        bwa_index        = ch_bwa_index
        bgzip_fasta_file = ch_bgzipped_fasta
        fasta_gzi_index  = ch_gzi_index
}
