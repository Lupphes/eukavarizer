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
include { DELLY_CALL                                } from '../modules/nf-core/delly/call/main.nf'
include { SAMPLE_REHEADER as DELLY_SAMPLE_REHEADER  } from '../modules/local/sample_regen/main.nf'
include { SVYNC as DELLY_SVYNC                      } from '../modules/nf-core/svync/main'
include { GUNZIP as DELLY_GUN                       } from '../modules/nf-core/gunzip/main'

/*
    ──────────────────────────────────────────────────────────
    NF-CORE MODULE IMPORTS (commented out for future use)
    ─────────────────────────────────────────────────────────
*/

/*
 * MANTA_GERMLINE
 *   For germline short-read SV calling (large deletions, inversions, etc.).
 */
include { MANTA_GERMLINE                                        } from '../modules/nf-core/manta/germline/main'
include { SAMPLE_REHEADER as MANTA_SMALL_SAMPLE_REHEADER        } from '../modules/local/sample_regen/main.nf'
include { SAMPLE_REHEADER as MANTA_CANDIDATE_SAMPLE_REHEADER    } from '../modules/local/sample_regen/main.nf'
include { SAMPLE_REHEADER as MANTA_DIPLOID_SAMPLE_REHEADER      } from '../modules/local/sample_regen/main.nf'
include { SVYNC as MANTA_SMALL_SVYNC                            } from '../modules/nf-core/svync/main'
include { SVYNC as MANTA_CANDIDATE_SVYNC                        } from '../modules/nf-core/svync/main'
include { SVYNC as MANTA_DIPLOID_SVYNC                          } from '../modules/nf-core/svync/main'

/*
 * GRIDSS_GRIDSS
 *   Advanced, local assembly-based SV caller for short-read data
 *   (both germline and somatic).
 */
include { GRIDSS_GRIDSS                             } from '../modules/nf-core/gridss/gridss/main'
include { SAMPLE_REHEADER as GRIDSS_SAMPLE_REHEADER } from '../modules/local/sample_regen/main.nf'
include { SVYNC as GRIDSS_SVYNC                     } from '../modules/nf-core/svync/main'
include { GUNZIP as GRIDSS_GUN                      } from '../modules/nf-core/gunzip/main'

/*
 * DYSGU
 *   Machine-learning approach for structural variant calling from short reads.
 */
include { DYSGU                                     } from '../modules/nf-core/dysgu/main'
include { SAMPLE_REHEADER as DYSGU_SAMPLE_REHEADER  }  from '../modules/local/sample_regen/main.nf'

/*
 * TIDDIT_SV
 *   Coverage-based short-read SV detection, used in some popular pipelines.
 */
include { TIDDIT_SV                                 } from '../modules/nf-core/tiddit/sv/main'
include { SAMPLE_REHEADER as TIDDIT_SAMPLE_REHEADER }  from '../modules/local/sample_regen/main.nf'

/*
 * SVABA
 *   Local assembly-based SV and indel detection (single-sample or tumor/normal).
 */
include { SVABA } from '../modules/nf-core/svaba/main' // Currently disabled

/*
 * SNIFFLES
 *   Long-read (PacBio/ONT) structural variant caller using split-read logic.
 */
include { SNIFFLES          } from '../modules/nf-core/sniffles/main' // Currently disabled
include { BCFTOOLS_FILTER   } from '../modules/nf-core/bcftools/filter/main'
include { BCFTOOLS_QUERY    } from '../modules/nf-core/bcftools/query/main'

/*
 * CUTESV
 *   Alternative long-read structural variant caller, similar scope to Sniffles.
 */
include { CUTESV                                    } from '../modules/nf-core/cutesv/main'
include { SAMPLE_REHEADER as CUTESV_SAMPLE_REHEADER }  from '../modules/local/sample_regen/main.nf'

// TO INCLUDE
// nf-core modules install deepvariant/rundeepvariant
// nf-core modules install duphold (LUMPY)
// https://github.com/etal/cnvkit

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
        delly_flag
        manta_flag
        gridss_flag
        dysgu_flag
        tiddit_flag
        svaba_flag
        sniffles_flag
        cutesv_flag
        debug_flag

    main:
        view("🔬 Running Structural Variant Calling Pipeline")

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        STEP 1) Prepare Inputs (For DELLY)
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        bam_inputs = ch_bam_files
            .join(ch_bam_indexes, by: 0)
            .map { meta, bam, bai -> tuple(meta, bam, bai) }

        fasta_input = ch_genome_file
            .map { fasta ->
                def id = fasta.baseName.replaceAll(/\.fna$/, '')
                tuple([id: id], fasta)
            }

        fai_input = ch_fasta_index
            .map { meta, fai -> tuple(meta, fai) }

        bwa_input = ch_bwa_index
            .map { idx -> tuple([id:'bwa_index'], idx) }

        if (debug_flag) {
            bam_inputs.view     { "DEBUG: BAM FILES -> ${it}" }
            fasta_input.view    { "DEBUG: FASTA -> ${it}" }
            fai_input.view      { "DEBUG: FAI -> ${it}" }
            bwa_input.view      { "DEBUG: BWA INDEX -> ${it}" }
        }

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        STEP 2) DELLY CALL
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        if (delly_flag) {
            name_delly = "delly"

            ch_delly_input = ch_bam_files
                .map { meta, bam ->
                    def meta_delly = meta + [id: "${meta.id}_${name_delly}", prefix: "${meta.prefix}_${name_delly}"]
                    tuple(meta_delly, bam)
                }
                .join(
                    ch_bam_indexes.map { meta, bai ->
                        tuple(meta + [id: "${meta.id}_${name_delly}", prefix: "${meta.prefix}_${name_delly}"], bai)
                    },
                    by: 0
                )
                .map { meta, bam, bai -> tuple(meta, bam, bai, [], [], []) }

            structural_delly = DELLY_CALL(ch_delly_input, fasta_input, fai_input)

            // Rename samples in the VCF
            structural_delly_reheaded = DELLY_SAMPLE_REHEADER(
                structural_delly.bcf.map { meta, vcf -> tuple(meta, vcf) },
                name_delly
            )

            // Change headers and generate GZ and UnGZFiles
            structural_delly_svync = DELLY_SVYNC(
                structural_delly_reheaded.vcf
                    .join(structural_delly_reheaded.tbi, by: 0)
                    .map { meta, vcf, tbi ->
                        tuple([id: "${meta.id}_svync"], vcf, tbi)
                    }
                    .combine(channel.value(file("${projectDir}/assets/svync/${name_delly}.yaml")))
            )

            // Unzip the SVYNC VCF
            structural_delly_svync_unzipped = DELLY_GUN(
                structural_delly_svync.vcf
            ).gunzip

            if (debug_flag) {
                structural_delly_reheaded.vcf.view      { "DEBUG: DELLY REHEADERED VCF -> ${it}" }
                structural_delly_reheaded.vcfgz.view    { "DEBUG: DELLY REHEADERED VCF.GZ -> ${it}" }
                structural_delly_reheaded.tbi.view      { "DEBUG: DELLY REHEADERED GZ TBI -> ${it}" }
                structural_delly_reheaded.csi.view      { "DEBUG: DELLY REHEADERED CSI -> ${it}" }
                structural_delly_svync_unzipped.view    { "DEBUG: DELLY SVYNC UNZIPPED VCF -> ${it}" }
                structural_delly_svync.vcf.view         { "DEBUG: DELLY SVYNC BGZIPPED VCF -> ${it}" }
                structural_delly_svync.tbi.view         { "DEBUG: DELLY SVYNC TBI INPUT -> ${it}" }
            }
        }

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        STEP 3) MANTA
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        if (manta_flag) {
            name_manta = "manta"

            structural_manta = MANTA_GERMLINE(
                bam_inputs.map { meta, bam, bai ->
                    tuple(meta + [id: "${meta.id}_${name_manta}", prefix: "${meta.prefix}_${name_manta}"], bam, bai, [], []) },
                fasta_input.map { meta, fasta ->
                    tuple(meta + [id: "${meta.id}_${name_manta}", prefix: "${meta.prefix}_${name_manta}"], fasta) },
                ch_fasta_index_gz.map { meta, fai ->
                    tuple(meta + [id: "${meta.id}_${name_manta}", prefix: "${meta.prefix}_${name_manta}"], fai) },
                []
            )

            // Rename samples in the VCF
            structural_manta_small_reheaded = MANTA_SMALL_SAMPLE_REHEADER(
                structural_manta.candidate_small_indels_vcf.map { meta, vcf -> tuple(meta, vcf) },
                "${name_manta}_small"
            )

            structural_manta_candidate_reheaded = MANTA_CANDIDATE_SAMPLE_REHEADER(
                structural_manta.candidate_sv_vcf.map { meta, vcf -> tuple(meta, vcf) },
                "${name_manta}_candidate"
            )

            structural_manta_diploid_reheaded = MANTA_DIPLOID_SAMPLE_REHEADER(
                structural_manta.diploid_sv_vcf.map { meta, vcf -> tuple(meta, vcf) },
                "${name_manta}_diploid"
            )

            if (debug_flag) {
                structural_manta_small_reheaded.vcf.view            { "DEBUG: MANTA SMALL REHEADERED VCF -> ${it}" }
                structural_manta_small_reheaded.vcfgz.view          { "DEBUG: MANTA SMALL REHEADERED VCF.GZ -> ${it}" }
                structural_manta_small_reheaded.tbi.view            { "DEBUG: MANTA SMALL REHEADERED GZ TBI -> ${it}" }
                structural_manta_small_reheaded.csi.view            { "DEBUG: MANTA SMALL REHEADERED CSI -> ${it}" }

                structural_manta_candidate_reheaded.vcf.view        { "DEBUG: MANTA CANDIDATE REHEADERED VCF -> ${it}" }
                structural_manta_candidate_reheaded.vcfgz.view      { "DEBUG: MANTA CANDIDATE REHEADERED VCF.GZ -> ${it}" }
                structural_manta_candidate_reheaded.tbi.view        { "DEBUG: MANTA CANDIDATE REHEADERED GZ TBI -> ${it}" }
                structural_manta_candidate_reheaded.csi.view        { "DEBUG: MANTA CANDIDATE REHEADERED CSI -> ${it}" }

                structural_manta_diploid_reheaded.vcf.view          { "DEBUG: MANTA DIPLOID REHEADERED VCF -> ${it}" }
                structural_manta_diploid_reheaded.vcfgz.view        { "DEBUG: MANTA DIPLOID REHEADERED VCF.GZ -> ${it}" }
                structural_manta_diploid_reheaded.tbi.view          { "DEBUG: MANTA DIPLOID REHEADERED GZ TBI -> ${it}" }
                structural_manta_diploid_reheaded.csi.view          { "DEBUG: MANTA DIPLOID REHEADERED CSI -> ${it}" }
            }
        }

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        STEP 4) GRIDSS_GRIDSS
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        if (gridss_flag) {
            name_gridss = "gridss"

            structural_gridss = GRIDSS_GRIDSS(
                bam_inputs.map { meta, bam, bai -> tuple(meta + [id: "${meta.id}_${name_gridss}", prefix: "${meta.prefix}_${name_gridss}"], bam) },
                ch_fasta_file.map { fasta -> tuple([id: fasta.baseName], fasta) },
                fai_input.map { meta, fai -> tuple(meta + [id: "${meta.id}_${name_gridss}", prefix: "${meta.prefix}_${name_gridss}"], fai) },
                bwa_input.map { meta, bwa -> tuple(meta + [id: "${meta.id}_${name_gridss}", prefix: "${meta.prefix}_${name_gridss}"], bwa[1]) }
            )

            // Rename samples in the VCF
            structural_gridss_reheaded = GRIDSS_SAMPLE_REHEADER(
                structural_gridss.vcf.map { meta, vcf -> tuple(meta, vcf) },
                name_gridss
            )

            // Change headers and generate GZ and UnGZFiles
            structural_gridss_svync = GRIDSS_SVYNC(
                structural_gridss_reheaded.vcf
                    .join(structural_gridss_reheaded.tbi, by: 0)
                    .map { meta, vcf, tbi ->
                        tuple([id: "${meta.id}_svync"], vcf, tbi)
                    }
                    .combine(channel.value(file("${projectDir}/assets/svync/${name_gridss}.yaml")))
            )

            // Unzip the SVYNC VCF
            structural_gridss_unzipped = GRIDSS_GUN(
                structural_gridss_svync.vcf
            ).gunzip

            if (debug_flag) {
                structural_gridss_reheaded.vcf.view     { "DEBUG: GRIDSS REHEADERED VCF -> ${it}" }
                structural_gridss_reheaded.vcfgz.view   { "DEBUG: GRIDSS REHEADERED VCF.GZ -> ${it}" }
                structural_gridss_reheaded.tbi.view     { "DEBUG: GRIDSS REHEADERED GZ TBI -> ${it}" }
                structural_gridss_reheaded.csi.view     { "DEBUG: GRIDSS REHEADERED CSI -> ${it}" }
                structural_gridss_unzipped.view         { "DEBUG: GRIDSS UNZIPPED VCF -> ${it}" }
                structural_gridss_svync.vcf.view        { "DEBUG: GRIDSS SVYNC BGZIPPED VCF -> ${it}" }
                structural_gridss_svync.tbi.view        { "DEBUG: GRIDSS SVYNC TBI INPUT -> ${it}" }
            }
        }

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        STEP 5) DYSGU
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        if (dysgu_flag) {
            name_dysgu = "dysgu"

            ch_dysgu_input_bam = ch_bam_files
                .map { meta, bam ->
                    def meta_dysgu = meta + [id: "${meta.id}_${name_dysgu}", prefix: "${meta.prefix}_${name_dysgu}"]
                    tuple(meta_dysgu, bam)
                }
                .join(
                    ch_bam_indexes.map { meta, bai ->
                        tuple(meta + [id: "${meta.id}_${name_dysgu}", prefix: "${meta.prefix}_${name_dysgu}"], bai)
                    },
                    by: 0
                )
                .map { meta, bam, bai -> tuple(meta, bam, bai) }

            ch_dysgu_input_fai = fasta_input
                .map { meta, fasta ->
                    def meta_dysgu = meta + [id: "${meta.id}_${name_dysgu}", prefix: "${meta.prefix}_${name_dysgu}"]
                    tuple(meta_dysgu, fasta)
                }
                .join(
                    fai_input.map { meta, fai ->
                        tuple(meta + [id: "${meta.id}_${name_dysgu}", prefix: "${meta.prefix}_${name_dysgu}"], fai)
                    },
                    by: 0
                )
                .map { meta, fasta, fai -> tuple(meta, fasta, fai) }



            structural_dysgu = DYSGU(
                ch_dysgu_input_bam,
                ch_dysgu_input_fai,
            )

            // Rename samples in the VCF
            structural_dysgu_reheaded = DYSGU_SAMPLE_REHEADER(
                structural_dysgu.vcf.map { meta, vcf -> tuple(meta, vcf) },
                name_dysgu
            )

            // Debug view for DYSGU
            if (debug_flag) {
                structural_dysgu_reheaded.vcf.view      { "DEBUG: DYSGU VCF INPUT -> ${it}" }
                structural_dysgu_reheaded.vcfgz.view    { "DEBUG: DYSGU BGZIPPED VCF -> ${it}" }
                structural_dysgu_reheaded.tbi.view      { "DEBUG: DYSGU TBI INPUT -> ${it}" }
                structural_dysgu_reheaded.csi.view      { "DEBUG: DYSGU CSI INPUT -> ${it}" }
            }
        }

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        STEP 6) TIDDIT_SV
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        if (tiddit_flag) {
            name_tiddit = "tiddit"

            structural_tiddit = TIDDIT_SV(
                bam_inputs.map { meta, bam, bai -> tuple(meta + [id: "${meta.id}_${name_tiddit}", prefix: "${meta.prefix}_${name_tiddit}"], bam, bai) },
                fasta_input,
                bwa_input.map { meta, bwa -> tuple(meta + [id: "${meta.id}_${name_tiddit}", prefix: "${meta.prefix}_${name_tiddit}"], bwa[1]) }
            )

            // Rename samples in the VCF
            structural_tiddit_reheaded = TIDDIT_SAMPLE_REHEADER(
                structural_tiddit.vcf.map { meta, vcf -> tuple(meta, vcf) },
                name_tiddit
            )

            if (debug_flag) {
                structural_tiddit_reheaded.vcf.view      { "DEBUG: TIDDIT VCF INPUT -> ${it}" }
                structural_tiddit.ploidy.view            { "DEBUG: TIDDIT PLOIDY INPUT -> ${it}" }
                structural_tiddit_reheaded.vcfgz.view    { "DEBUG: TIDDIT BGZIPPED VCF -> ${it}" }
                structural_tiddit_reheaded.tbi.view      { "DEBUG: TIDDIT TBI INPUT -> ${it}" }
                structural_tiddit_reheaded.csi.view      { "DEBUG: TIDDIT CSI INPUT -> ${it}" }
            }
        }

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        STEP 7) SVABA
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        if (svaba_flag) {
            // structural_svaba = SVABA(
            //     ch_bam_files
            //         .join(ch_bam_indexes, by: 0)
            //         .map { meta, bam, bai -> tuple(meta, null, null, bam, bai) },
            //     ch_genome_file.map { fasta -> tuple([id: fasta.baseName], fasta) },
            //     ch_fasta_index.map { meta, fai -> tuple(meta, fai) },
            //     ch_bwa_index.map { meta, bwa -> tuple(meta, bwa) },
            //     "",
            //     "",
            //     ""
            // )
            if (debug_flag) {
                // structural_svaba.view { "DEBUG: SVABA INPUT -> ${it}" }
            }
        }

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        STEP 8) SNIFFLES
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        if (sniffles_flag) {
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
            if (debug_flag) {
                // structural_sniffles.view { "DEBUG: SNIFFLES INPUT -> ${it}" }
            }
        }

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        STEP 9) CUTESV
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        if (cutesv_flag) {
            name_cutesv = "cutesv"

            structural_cutesv = CUTESV(
                bam_inputs.map { meta, bam, bai ->
                    tuple(meta + [id: "${meta.id}_${name_cutesv}"], bam, bai) },
                fasta_input.map { meta, fasta -> tuple(meta + [id: "${meta.id}_${name_cutesv}"], fasta) }
            )

            // Rename samples in the VCF
            structural_cutesv_reheaded = CUTESV_SAMPLE_REHEADER(
                structural_cutesv.vcf.map { meta, vcf -> tuple(meta, vcf) },
                name_cutesv
            )

            if (debug_flag) {
                structural_cutesv_reheaded.vcf.view      { "DEBUG: CUTESV VCF INPUT -> ${it}" }
                structural_cutesv_reheaded.vcfgz.view    { "DEBUG: CUTESV BGZIPPED VCF -> ${it}" }
                structural_cutesv_reheaded.tbi.view      { "DEBUG: CUTESV TBI INPUT -> ${it}" }
                structural_cutesv_reheaded.csi.view      { "DEBUG: CUTESV CSI INPUT -> ${it}" }
            }
        }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        STEP 10) Emit Outputs
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    emit:
        delly_variants                      = delly_flag  ? (structural_delly_svync_unzipped                ?: null) : null
        delly_variants_index                = delly_flag  ? (structural_delly_svync?.tbi                    ?: null) : null
        delly_variants_bgzipped             = delly_flag  ? (structural_delly_svync?.vcf                    ?: null) : null
        delly_variants_index_bgzipped       = manta_flag  ? (structural_delly_reheaded?.vcfgz               ?: null) : null

        manta_small_variants                = manta_flag  ? (structural_manta_small_reheaded?.vcf           ?: null) : null
        manta_small_variants_index          = manta_flag  ? (structural_manta_small_reheaded?.tbi           ?: null) : null
        manta_small_variants_bgzipped       = manta_flag  ? (structural_manta_small_reheaded?.vcfgz         ?: null) : null

        manta_candidate_variants            = manta_flag  ? (structural_manta_candidate_reheaded?.vcf       ?: null) : null
        manta_candidate_variants_index      = manta_flag  ? (structural_manta_candidate_reheaded?.tbi       ?: null) : null
        manta_candidate_variants_bgzipped   = manta_flag ? (structural_manta_candidate_reheaded?.vcfgz     ?: null) : null

        manta_diploid_variants              = manta_flag  ? (structural_manta_diploid_reheaded?.vcf         ?: null) : null
        manta_diploid_variants_index        = manta_flag  ? (structural_manta_diploid_reheaded?.tbi         ?: null) : null
        manta_diploid_variants_bgzipped     = manta_flag  ? (structural_manta_diploid_reheaded?.vcfgz       ?: null) : null

        gridss_variants                     = gridss_flag ? (structural_gridss_unzipped                     ?: null) : null
        gridss_variants_index               = gridss_flag ? (structural_gridss_svync?.tbi                   ?: null) : null
        gridss_variants_bgzipped            = gridss_flag ? (structural_gridss_svync?.vcf                   ?: null) : null

        dysgu_variants                      = dysgu_flag  ? (structural_dysgu_reheaded?.vcf                 ?: null) : null
        dysgu_variants_index                = dysgu_flag  ? (structural_dysgu_reheaded?.tbi                 ?: null) : null
        dysgu_variants_bgzipped             = dysgu_flag  ? (structural_dysgu_reheaded?.vcfgz               ?: null) : null

        tiddit_variants                     = tiddit_flag ? (structural_tiddit_reheaded?.vcf                ?: null) : null
        tiddit_ploidy                       = tiddit_flag ? (structural_tiddit?.ploidy                      ?: null) : null
        tiddit_variants_index               = tiddit_flag ? (structural_tiddit_reheaded?.tbi                ?: null) : null
        tiddit_variants_bgzipped            = tiddit_flag ? (structural_tiddit_reheaded?.vcfgz              ?: null) : null

        cutesv_variants                     = cutesv_flag ? (structural_cutesv_reheaded?.vcf                ?: null) : null
        cutesv_variants_index               = cutesv_flag ? (structural_cutesv_reheaded?.tbi                ?: null) : null
        cutesv_variants_bgzipped            = cutesv_flag ? (structural_cutesv_reheaded?.vcfgz              ?: null) : null
}
