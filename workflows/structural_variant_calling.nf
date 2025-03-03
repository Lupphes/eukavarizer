/*
 * Pipeline: STRUCTURAL_VARIANT_CALLING
 * This workflow demonstrates a minimal SV calling pipeline
 * that currently runs only DELLY. Future expansions can enable
 * Manta, GRIDSS, Dysgu, TIDDIT, SVABA, Sniffles, and CuteSV.
 */

/*
    ──────────────────────────────────────────────────────────
    LOCAL MODULE IMPORTS
    ──────────────────────────────────────────────────────────
*/

/*
 * DELLY_CALL
 *   Local module that calls structural variants (deletions, duplications,
 *   inversions, translocations) in short-read data using Delly.
 */
include { DELLY_CALL } from '../modules/local/delly/call/main.nf'


/*
    ──────────────────────────────────────────────────────────
    NF-CORE MODULE IMPORTS (commented out for future use)
    ─────────────────────────────────────────────────────────
*/

/*
 * MANTA_GERMLINE
 *   For germline short-read SV calling (large deletions, inversions, etc.).
 */
include { MANTA_GERMLINE } from '../modules/nf-core/manta/germline/main'

/*
 * GRIDSS_GRIDSS
 *   Advanced, local assembly-based SV caller for short-read data
 *   (both germline and somatic).
 */
// include { GRIDSS_GRIDSS } from '../modules/nf-core/gridss/gridss/main'

/*
 * DYSGU
 *   Machine-learning approach for structural variant calling from short reads.
 */
// include { DYSGU } from '../modules/nf-core/dysgu/main'

/*
 * TIDDIT_SV
 *   Coverage-based short-read SV detection, used in some popular pipelines.
 */
// include { TIDDIT_SV } from '../modules/nf-core/tiddit/sv/main'

/*
 * SVABA
 *   Local assembly-based SV and indel detection (single-sample or tumor/normal).
 */
// include { SVABA } from '../modules/nf-core/svaba/main'

/*
 * SNIFFLES
 *   Long-read (PacBio/ONT) structural variant caller using split-read logic.
 */
// include { SNIFFLES } from '../modules/nf-core/sniffles/main'

/*
 * CUTESV
 *   Alternative long-read structural variant caller, similar scope to Sniffles.
 */
// include { CUTESV } from '../modules/nf-core/cutesv/main'


/*
    ──────────────────────────────────────────────────────────
    WORKFLOW DEFINITION
    ──────────────────────────────────────────────────────────
*/

workflow STRUCTURAL_VARIANT_CALLING {

    /*
     * Input channels for:
     *   1) ch_bam_files    => Aligned BAM files
     *   2) ch_bam_indexes  => Corresponding BAM/CRAM indices
     *   3) ch_genome_file  => Reference genome (FASTA)
     *   4) ch_fasta_index  => FASTA index (.fai)
     *   5) ch_bwa_index    => Directory containing BWA index files
     */
    take:
        ch_bam_files
        ch_bam_indexes
        ch_genome_file
        ch_fasta_index
        ch_bwa_index

    main:
        view("🔬 Running Structural Variant Calling (currently only DELLY)")

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        STEP 1) Prepare Inputs (For DELLY)
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        // A) Bam + BAI -> (meta, bam, bai)
        ch_bam_files.view   { "DEBUG: BAM FILES -> ${it}" }
        ch_bam_indexes.view { "DEBUG: BAI FILES -> ${it}" }
        bam_inputs = ch_bam_files
            .join(ch_bam_indexes, by: 0)
            .map { meta, bam, bai -> tuple(meta, bam, bai) }
            .view { "DEBUG: DELLY input -> ${it}" }

        // B) FASTA -> (meta2, fasta)
        ch_genome_file.view { "DEBUG: FASTA -> ${it}" }
        fasta_input = ch_genome_file
            .map { fasta -> tuple([id: fasta.baseName], fasta) }
            .view { "DEBUG: DELLY FASTA -> ${it}" }

        // C) FAI -> (meta3, fai)
        ch_fasta_index.view { "DEBUG: FAI -> ${it}" }
        fai_input = ch_fasta_index
            .map { meta, fai -> tuple(meta, fai) }
            .view { "DEBUG: DELLY FAI -> ${it}" }

        // D) BWA index -> (meta4, directory) - Not used by DELLY, but
        // may be required by GRIDSS or TIDDIT in the future
        ch_bwa_index.view { "DEBUG: BWA INDEX -> ${it}" }
        bwa_input = ch_bwa_index
            .map { idx -> tuple([id:'bwa_index'], idx) }
            .view { "DEBUG: BWA input -> ${it}" }

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        STEP 2) DELLY CALL (Active)
        - Calls short-read SVs (deletions, duplications, inversions, etc.)
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        structural_delly = DELLY_CALL(
            bam_inputs,     // (meta, bam, bai)
            fasta_input,    // (meta2, fasta)
            fai_input       // (meta3, fai)
        )

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        STEP 3) MANTA_GERMLINE (Commented)
        - For germline short-read SVs, including large deletions, inversions, etc.
        - Typically: (meta,bam,bai, [bed?],[bed?]), (meta2,fasta), (meta3,fai)
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        // structural_manta_germline = MANTA_GERMLINE(
        //     bam_inputs,
        //     fasta_input,
        //     fai_input,
        //     null  // or Channel.empty() if you have no intervals
        // )


        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        STEP 4) GRIDSS_GRIDSS (Commented)
        - Advanced local assembly-based approach for short-read data,
        capturing complex rearrangements (germline or somatic).
        - Usually: (meta,bam), (meta2,fasta), (meta3,fai), (meta4,bwa_index).
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        // structural_gridss = GRIDSS_GRIDSS(
        //     bam_inputs.map { meta, bam, bai -> tuple(meta, bam) },
        //     fasta_input,
        //     fai_input,
        //     bwa_input
        // )


        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        STEP 5) DYSGU (Commented)
        - A machine-learning approach for short-read SV detection.
        - Typically: (meta,bam,bai), (meta2,fasta), (meta3,fai).
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        // structural_dysgu = DYSGU(
        //     bam_inputs,
        //     fasta_input,
        //     fai_input
        // )


        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        STEP 6) TIDDIT_SV (Commented)
        - Coverage-based short-read SV detection method used in some pipelines.
        - Typically: (meta,bam,bai), (meta2,fasta), (meta3,bwa_index?).
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        // structural_tiddit = TIDDIT_SV(
        //     bam_inputs,
        //     fasta_input,
        //     bwa_input // if TIDDIT needs BWA index
        // )


        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        STEP 7) SVABA (Commented)
        - Local assembly-based method to detect SVs & Indels in single-sample
        or tumor/normal short-read data.
        - Minimal usage: (meta,bam,bai), (meta2,fasta).
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        // structural_svaba = SVABA(
        //     bam_inputs,
        //     fasta_input
        // )


        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        STEP 8) SNIFFLES (Commented)
        - For long-read (PacBio/ONT) SV calling using split-read signals.
        - Typically 5 arguments: (meta,bam,bai), (meta2,fasta),
        (meta3,tandem_file?), val vcf_output, val snf_output.
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        // structural_sniffles = SNIFFLES(
        //     bam_inputs,
        //     fasta_input,
        //     Channel.value( tuple([id:'dummy'], file('tandem.bed')) ),
        //     Channel.value("sniffles_output.vcf.gz"),
        //     Channel.value("sniffles_output.snf")
        // )


        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        STEP 9) CUTESV (Commented)
        - Another long-read SV caller, similar in scope to Sniffles.
        - Typically: (meta,bam,bai), (meta2,fasta), (meta3,fai).
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        // structural_cutesv = CUTESV(
        //     bam_inputs,
        //     fasta_input,
        //     fai_input
        // )



    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        STEP 10) Emit Outputs
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    emit:
        // Only DELLY is currently active -> BCF + CSI index
        delly_variants       = structural_delly.bcf
        delly_variants_index = structural_delly.csi

    // For future expansions, you'd add Manta, GRIDSS, etc. outputs below:
    // manta_germline_sv_vcf = structural_manta_germline.candidate_sv_vcf
    // ...
}
