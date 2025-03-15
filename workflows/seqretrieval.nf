/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { REFSEQ_RETRIEVER      } from '../modules/local/refseq_retriever/main'
include { SEQUENCE_RETRIEVER    } from '../modules/local/sequence_retriever/main'

// nf-core modules install pbtk/bam2fastq

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
        // STEP 1: Retrieve Reference Genome (Outputs JSON)
        //
        ch_refseq_json = REFSEQ_RETRIEVER(
            ch_taxonomy_id,
            ch_outdir,
            ch_genome_file
        )

        parsed_ch_refseq_json = ch_refseq_json.refseq_json
            .splitJson()
            .collect()
            .map { jsonList ->
                def jsonMap = jsonList.collectEntries { [it.key, it.value] }
                return tuple(jsonMap['genome_size'] as Integer, jsonMap['genome_size_ungapped'] as Integer, file(jsonMap['genome_file']))
            }

        //
        // STEP 3: Retrieve Sequencing Data (Outputs JSON)
        //
        ch_ena_results = SEQUENCE_RETRIEVER(
            ch_taxonomy_id,
            parsed_ch_refseq_json.map { it[1] },
            ch_outdir,
            ch_library_strategy.map { it.join(' ').trim() },
            ch_instrument_platform.map { it.join(' ').trim() },
            ch_minimum_coverage,
            ch_maximum_coverage,
            ch_max_results,
            ch_assembly_quality,
            ch_sequence_dir
        )

        raw_fastqs = ch_ena_results.fastq_files

        if (debug_flag) {
            raw_fastqs.view { "DEBUG: raw_fastqs -> ${it.getClass()} | ${it}" }
        }

        // Correct file grouping with proper handling of ArrayList
        grouped_fastqs = raw_fastqs
            .flatMap { file -> // Flatten in case of nested lists
                file instanceof List ? file : [file] // Ensure we are iterating over individual file paths
            }
            .map { file -> tuple(file.getParent().getName(), file) } // Extract run_accession from parent folder
            .groupTuple()
            .map { run_accession, files ->
                def paired = files.findAll { it.toString().contains('_1') || it.toString().contains('_2') }
                def unpaired = files - paired
                tuple(run_accession, paired, unpaired)
            }

        if (debug_flag) {
            grouped_fastqs.view { "DEBUG: grouped_fastqs -> ${it}" }
        }

    emit:
        genome_file = ch_refseq_json.genome_file // genome_file
        genome_size_ungapped = parsed_ch_refseq_json.map { it[1] } // genome_size_ungapped
        fastq_files = ch_ena_results.fastq_files // fullPaths
        grouped_fastqs = grouped_fastqs // grouped_fastqs
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
