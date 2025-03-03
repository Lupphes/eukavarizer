include { DELLY_CALL } from '../modules/local/delly/call/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW: STRUCTURAL_VARIANT_CALLING
    - Identifies structural variants using DELLY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow STRUCTURAL_VARIANT_CALLING {

    take:
        ch_bam_files     // Aligned BAM files
        ch_bam_indexes   // BAM index files (.bai)
        ch_genome_file   // Reference genome (FASTA)
        ch_fasta_index   // FASTA index (.fai)

    main:
        log.info "🔬 Running Structural Variant Calling with DELLY"

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            STEP 1: Prepare BAM Input for DELLY
            - Joins BAM files with their respective index (.bai)
            - Converts to a 3-element tuple: (meta, bam, bai)
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        def delly_input = ch_bam_files
            .join(ch_bam_indexes, by: 0)
            .map { meta, bam, bai -> tuple(meta, bam, bai) }
            .view { "DEBUG: DELLY BAM input -> ${it}" }

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            STEP 2: Prepare Reference Genome (FASTA)
            - Converts FASTA file into a 2-element tuple: (meta, fasta)
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        def delly_fasta = ch_genome_file
            .map { fasta -> tuple([id: fasta.baseName], fasta) }
            .view { "DEBUG: DELLY FASTA input -> ${it}" }

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            STEP 3: Prepare FASTA Index (.fai)
            - Converts FASTA index into a 2-element tuple: (meta, fai)
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        def delly_fai = ch_fasta_index
            .map { meta, fai -> tuple(meta, fai) }
            .view { "DEBUG: DELLY FASTA index -> ${it}" }

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            STEP 4: Run DELLY for Structural Variant Calling
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        structural_results = DELLY_CALL(
            delly_input,   // BAM + index tuple
            delly_fasta,   // FASTA tuple
            delly_fai      // FASTA index tuple
        )

    emit:
        variants       = structural_results.bcf   // Output BCF variant file
        variants_index = structural_results.csi   // Output index file for BCF
}
