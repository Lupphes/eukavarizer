/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { REFSEQ_RETRIEVER       } from '../modules/local/refseq_retriver/main'
include { SEQUENCE_RETRIEVER     } from '../modules/local/sequence_retriver/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SEQRETRIEVAL {

    take:
        ch_taxonomy_id
        ch_outdir
        ch_local_sequences_dir
        ch_local_refseq_path

    main:
        //
        // STEP 1: Retrieve Reference Genome (Outputs JSON)
        //
        ch_refseq_json = REFSEQ_RETRIEVER(
            ch_taxonomy_id,
            ch_outdir,
            ch_local_refseq_path
        ).out.refseq_json

        //
        // STEP 2: Parse Genome JSON
        //
        parsed_refseq_output = ch_refseq_json
            .splitJson()
            .map { json -> tuple(json.genome_size, json.genome_size_ungapped, file(json.genome_file)) }

        //
        // STEP 3: Retrieve Sequencing Data (Outputs JSON)
        //
        ch_sequence_json = SEQUENCE_RETRIEVER(
            ch_taxonomy_id,
            parsed_refseq_output.map { it[1] }, // genome_size_ungapped
            ch_outdir
        ).out.sequence_json

        //
        // STEP 4: Parse Sequencing Data JSON
        //
        parsed_sequence_output = ch_sequence_json
            .splitJson()
            .map { json -> json.fastq_files.collect { file(it) } }

    emit:
        refseq_path = parsed_refseq_output.map { it[2] } // genome_file
        genome_size_ungapped = parsed_refseq_output.map { it[1] }
        fastq_files = parsed_sequence_output
}
