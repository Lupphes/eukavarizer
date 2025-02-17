/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { REFSEQ_RETRIEVER      } from '../modules/local/refseq_retriever/main'
include { SEQUENCE_RETRIEVER    } from '../modules/local/sequence_retriever/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SEQRETRIEVAL {

    take:
        ch_taxonomy_id
        ch_outdir
        ch_sequences_dir
        ch_local_refseq_path

    main:
        ch_library_strategy    = params.library_strategy ? Channel.value(params.library_strategy) : Channel.value("")
        ch_instrument_platform = params.instrument_platform ? Channel.value(params.instrument_platform) : Channel.value("")
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
            ch_local_refseq_path
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
            ch_library_strategy,
            ch_instrument_platform,
            ch_minimum_coverage,
            ch_maximum_coverage,
            ch_max_results,
            ch_assembly_quality,
            ch_sequences_dir
        )

    emit:
        refseq_path = ch_refseq_json.genome_file // genome_file
        genome_size_ungapped = parsed_ch_refseq_json.map { it[1] } // genome_size_ungapped
        fastq_files = ch_ena_results.sequences // fullPaths
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
