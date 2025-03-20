/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BIODBCORE_REFSEQ      } from '../modules/local/refseq_retriever/main'
include { SEQUENCE_RETRIEVER    } from '../modules/local/sequence_retriever/main'

include { SAMTOOLS_COLLATEFASTQ as BAM_SAMTOOLS_COLLATEFASTQ } from '../modules/nf-core/samtools/collatefastq/main'
include { SAMTOOLS_COLLATEFASTQ as CRAM_SAMTOOLS_COLLATEFASTQ } from '../modules/nf-core/samtools/collatefastq/main'

include { GUNZIP         }      from '../modules/nf-core/gunzip/main'
include { TABIX_BGZIP    }      from '../modules/nf-core/tabix/bgzip/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SEQRETRIEVAL {

    take:
        ch_taxonomy_id
        ch_outdir
        ch_sequence_dir
        ch_genome_file
        debug_flag

    main:
        ch_library_strategy    = params.library_strategy ? Channel.value(params.library_strategy) : Channel.value([])
        ch_instrument_platform = params.instrument_platform ? Channel.value(params.instrument_platform) : Channel.value([])
        ch_minimum_coverage    = params.minimum_coverage ? Channel.value(params.minimum_coverage) : Channel.value(1)
        ch_maximum_coverage    = params.maximum_coverage ? Channel.value(params.maximum_coverage) : Channel.value("")
        ch_max_results         = params.max_results ? Channel.value(params.max_results) : Channel.value(1)
        ch_assembly_quality    = params.assembly_quality ? Channel.value(params.assembly_quality) : Channel.value("")

        //
        // STEP 1: Retrieve Reference Genome
        //
        refseq = BIODBCORE_REFSEQ(
            ch_taxonomy_id,
            ch_outdir,
            ch_genome_file
        )

        parsed_ch_refseq_json = refseq.json
            .splitJson()
            .collect()
            .map { jsonList ->
                def jsonMap = jsonList.collectEntries { [it.key, it.value] }
                return tuple(jsonMap['genome_size'] as Integer, jsonMap['genome_size_ungapped'] as Integer, file(jsonMap['genome_file']))
            }

        //
        // STEP 3: Retrieve Sequencing Data
        //
        ena = SEQUENCE_RETRIEVER(
            ch_taxonomy_id,
            parsed_ch_refseq_json.map   { it[1] },
            ch_outdir,
            ch_library_strategy.map     { it.join(' ').trim() },
            ch_instrument_platform.map  { it.join(' ').trim() },
            ch_minimum_coverage,
            ch_maximum_coverage,
            ch_max_results,
            ch_assembly_quality,
            ch_sequence_dir
        )

        raw_fastqs = ena.fastq_files
        raw_bam = ena.bam_files
        raw_cram = ena.cram_files

        // raw_bam.view { "DEBUG: raw_bam -> ${it}" }
        // raw_cram.view { "DEBUG: raw_cram -> ${it}" }

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
                def unpaired = files - paired
                tuple(run_accession, paired, unpaired)
            }

        grouped_bam = raw_bam
            .flatMap { file -> file instanceof List ? file : [file] }
            .map { file -> tuple([id: file.getParent().getName()], file.toString()) }
            .groupTuple()

        grouped_cram = raw_cram
            .flatMap { file -> file instanceof List ? file : [file] }
            .map { file -> tuple([id: file.getParent().getName()], file.toString()) }
            .groupTuple()

        // grouped_fastqs.view { "DEBUG: grouped_fastqs -> ${it}" }
        // grouped_bam.view    { "DEBUG: grouped_bam -> ${it}" }
        // grouped_cram.view   { "DEBUG: grouped_cram -> ${it}" }



        // repete, remove //
        genome_unzipped = GUNZIP(
            refseq.reference_genome.map { file -> tuple([id: file.simpleName.replaceFirst(/\.gz$/, '')], file) }
        ).gunzip


        ch_bgzipped_fasta_tuple = TABIX_BGZIP(
            genome_unzipped.map { meta, fasta -> tuple(meta, fasta) }
        )

        ch_bgzipped_fasta = ch_bgzipped_fasta_tuple.output.map { meta, file -> file }
        ch_gzi_index = ch_bgzipped_fasta_tuple.gzi



        ref_meta = ch_bgzipped_fasta.map { file -> tuple([id: file.simpleName.replaceFirst(/\.gz$/, '')], file) }
        // ref_meta.view { "DEBUG: refseq.reference_genome -> ${it}" }

 // repete, remove //

        bam =  BAM_SAMTOOLS_COLLATEFASTQ(
            grouped_bam,
            ref_meta,
            false
        )

        cram = CRAM_SAMTOOLS_COLLATEFASTQ(
            grouped_cram,
            ref_meta,
            false
        )

        // Process BAM reads (paired + unpaired)
        bam_paired = bam.fastq.map { meta, fastqs ->
            tuple("${meta.id}_paired", fastqs, [])
        }


        bam_unpaired = bam.fastq_singleton.map { meta, fastqs ->
            tuple("${meta.id}_unpaired", [], fastqs ? [fastqs].flatten() : [])
        }

        // Process CRAM reads (paired + unpaired)
        cram_paired = cram.fastq.map { meta, fastqs ->
            tuple("${meta.id}_paired", fastqs, [])
        }

        cram_unpaired = cram.fastq_singleton.map { meta, fastqs ->
            tuple("${meta.id}_unpaired", [], fastqs ? [fastqs].flatten() : [])
        }

        // bam_unpaired.view { "DEBUG: bam_unpaired -> ${it}" }
        // cram_unpaired.view { "DEBUG: cram_unpaired -> ${it}" }



        // Merge into existing grouped_fastqs
        grouped_fastqs = grouped_fastqs
            .mix(bam_paired)
            .mix(bam_unpaired)
            .mix(cram_paired)
            .mix(cram_unpaired)

        // Debug to confirm consistent output
        // grouped_fastqs.view { "DEBUG: grouped_fastqs -> ${it}" }


        if (debug_flag) {
            raw_fastqs.view { "DEBUG: raw_fastqs -> ${it.getClass()} | ${it}" }
        }



    emit:
        genome_file                 = refseq.reference_genome              // genome_file
        genome_size_ungapped        = parsed_ch_refseq_json.map { it[1] }  // genome_size_ungapped
        fastq_files                 = ena.fastq_files                      // fullPaths
        grouped_fastqs              = grouped_fastqs                       // grouped_fastqs
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
