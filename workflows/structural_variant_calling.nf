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
include { DELLY_CALL } from '../modules/nf-core/delly/call/main.nf'

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
include { GRIDSS_GRIDSS } from '../modules/nf-core/gridss/gridss/main'

/*
 * DYSGU
 *   Machine-learning approach for structural variant calling from short reads.
 */
include { DYSGU } from '../modules/nf-core/dysgu/main'

/*
 * TIDDIT_SV
 *   Coverage-based short-read SV detection, used in some popular pipelines.
 */
include { TIDDIT_SV } from '../modules/nf-core/tiddit/sv/main'

/*
 * SVABA
 *   Local assembly-based SV and indel detection (single-sample or tumor/normal).
 */
include { SVABA } from '../modules/nf-core/svaba/main' // Currently disabled

/*
 * SNIFFLES
 *   Long-read (PacBio/ONT) structural variant caller using split-read logic.
 */
include { SNIFFLES } from '../modules/nf-core/sniffles/main' // Currently disabled
include { BCFTOOLS_FILTER } from '../modules/nf-core/bcftools/filter/main'
include { BCFTOOLS_QUERY } from '../modules/nf-core/bcftools/query/main'

/*
 * CUTESV
 *   Alternative long-read structural variant caller, similar scope to Sniffles.
 */
include { CUTESV } from '../modules/nf-core/cutesv/main'

/*
    ──────────────────────────────────────────────────────────
    WORKFLOW DEFINITION
    ──────────────────────────────────────────────────────────
*/

workflow STRUCTURAL_VARIANT_CALLING {

    /*
     * Input channels for:
     *   1) ch_bam_files        => Aligned BAM files
     *   2) ch_bam_indexes      => Corresponding BAM/CRAM indices
     *   3) ch_genome_file      => Reference genome (FASTA)
     *   4) ch_fasta_index      => FASTA index (.fai)
     *   5) ch_bwa_index        => Directory containing BWA index files
     *   6) ch_fasta_index_gz   => FASTA index (.fai) for gzipped FASTA
     *   7) ch_fasta_file       => FASTA file (uncompressed)
     */
    take:
        ch_bam_files
        ch_bam_indexes
        ch_genome_file
        ch_fasta_index
        ch_bwa_index
        ch_fasta_index_gz
        ch_fasta_file

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
        STEP 2) DELLY CALL
        - Calls short-read SVs (deletions, duplications, inversions, etc.)
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        structural_delly = DELLY_CALL(
            ch_bam_files
            .join(ch_bam_indexes, by: 0)
            .map { meta, bam, bai -> tuple(meta, bam, bai, [], [], []) },
            fasta_input,
            fai_input
        )

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        STEP 3) MANTA_GERMLINE
        - For germline short-read SVs, including large deletions, inversions, etc.
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        structural_manta_germline = MANTA_GERMLINE(
            bam_inputs.map { meta, bam, bai -> tuple(meta, bam, bai, [], []) },
            fasta_input,
            ch_fasta_index_gz,
            []
        )

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        STEP 4) GRIDSS_GRIDSS
        - Advanced local assembly-based approach for short-read data,
        capturing complex rearrangements (germline or somatic).
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        structural_gridss = GRIDSS_GRIDSS(
            bam_inputs.map { meta, bam, bai -> tuple(meta, bam) },
            ch_fasta_file.map { fasta -> tuple([id: fasta.baseName], fasta) },
            fai_input.map { meta, fai -> tuple(meta, fai) },
            bwa_input.map { meta, bwa -> tuple(meta, bwa[1]) }
        )

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        STEP 5) DYSGU
        - A machine-learning approach for short-read SV detection.
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */


        // ch_bam_files.view { it -> println "1 BAM FILE INPUT: $it" }
        // ch_bam_indexes.view { it -> println "1 BAM INDEX INPUT: $it" }
        // fasta_input.view { it -> println "1 FASTA INPUT: $it" }
        // ch_fasta_index.view { it -> println "1 FASTA INDEX INPUT: $it" }

        // structural_dysgu = DYSGU(
        //     ch_bam_files
        //         .join(ch_bam_indexes, by: 0)
        //         .map { meta, bam, bai -> tuple(meta, bam, bai) },
        //     fasta_input
        //         .map { meta, fasta -> tuple(fasta.simpleName, meta, fasta) } // Ensure consistent FASTA ID
        //         .join(
        //             ch_fasta_index.map { meta, fai -> tuple(fai.simpleName, meta, fai) }, // Ensure consistent FAI ID
        //             by: 0
        //         )
        //         .map { id, meta_fasta, fasta, meta_fai, fai -> tuple(meta_fasta, fasta, fai) } // Reassemble metadata
        // )



        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        STEP 6) TIDDIT_SV (Commented)
        - Coverage-based short-read SV detection method used in some pipelines.
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        structural_tiddit = TIDDIT_SV(
            bam_inputs,
            fasta_input,
            bwa_input.map { meta, bwa -> tuple(meta, bwa[1]) }
        )
        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        STEP 7) SVABA
        - Local assembly-based method to detect SVs & Indels in single-sample
        or tumor/normal short-read data.
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        // 🚀 Run SVABA
        // structural_svaba = SVABA(
        //     ch_bam_files
        //         .join(ch_bam_indexes, by: 0)
        //         .map { meta, bam, bai -> tuple(meta, null, null, bam, bai) },       // (meta, tumorbam, tumorbai, normalbam=null, normalbai=null)
        //     ch_genome_file.map { fasta -> tuple([id: fasta.baseName], fasta) },     // (meta2, fasta)
        //     ch_fasta_index.map { meta, fai -> tuple(meta, fai) },                   // (meta2, fasta_fai)
        //     ch_bwa_index.map { meta, bwa -> tuple(meta, bwa) },                     // (meta3, bwa_index)
        //     "",                                                                      // (meta4, dbsnp)
        //     "",                                                        // (meta4, dbsnp_tbi)
        //     ""                                                    // (meta5, regions)
        // )


        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        STEP 8) SNIFFLES
        - For long-read (PacBio/ONT) SV calling using split-read signals.
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        // filtered_tandem_vcf = BCFTOOLS_FILTER(
        //     structural_manta_germline.candidate_sv_vcf.map { meta, vcf, tbi -> tuple(meta, vcf, tbi) }
        // )

        // tandem_repeats_bed = BCFTOOLS_QUERY(
        //     filtered_tandem_vcf.vcf
        //     .join(filtered_tandem_vcf.tbi, by: 0)
        //     .map { meta, vcf, tbi -> tuple(meta, vcf, tbi) },
        //     [],
        //     [],
        //     []
        // )

        // structural_sniffles = SNIFFLES(
        //     ch_bam_files
        //         .join(ch_bam_indexes, by: 0)
        //         .map { meta, bam, bai -> tuple(meta, bam, bai) },
        //     ch_genome_file
        //         .map { fasta -> tuple([id: fasta.baseName], fasta) },
        //     tandem_repeats_bed,
        //     [],
        //     []
        // )

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        STEP 9) CUTESV (Commented)
        - Another long-read SV caller, similar in scope to Sniffles.
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        structural_cutesv = CUTESV(
            bam_inputs,
            fasta_input
        )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        STEP 10) Emit Outputs
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    emit:
        delly_variants                      = structural_delly.bcf
        delly_variants_index                = structural_delly.csi
        manta_sv_variants                   = structural_manta_germline.candidate_sv_vcf
        manta_sv_variants_index             = structural_manta_germline.candidate_sv_vcf_tbi
        manta_diploid_sv_variants           = structural_manta_germline.diploid_sv_vcf
        manta_diploid_sv_variants_index     = structural_manta_germline.diploid_sv_vcf_tbi
        gridss_variants                     = structural_gridss.vcf
        // dysgu_variants                      = structural_dysgu.vcf
        // dysgu_tbi                           = structural_dysgu.tbi
        tiddit_variants                     = structural_tiddit.vcf
        tiddit_ploidy                       = structural_tiddit.ploidy
        cutesv_variants                     = structural_cutesv.vcf
}
