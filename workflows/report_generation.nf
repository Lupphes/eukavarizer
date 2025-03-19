/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES FOR STRUCTURAL VARIANT MERGING AND SYNCHRONIZATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SURVIVOR_MERGE        }       from '../modules/nf-core/survivor/merge/main'
include { SURVIVOR_STATS        }       from '../modules/nf-core/survivor/stats/main'
include { BCFTOOLS_MERGE        }       from '../modules/nf-core/bcftools/merge/main'
include { VCF_REPGEN            }       from '../modules/local/vcf_repgen/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow REPORT_GENERATION {

    // Input: Structural variant files from various tools
    take:
        genome_file
        delly_variants
        delly_variants_bgzipped
        delly_variants_index
        delly_variants_index_bgzipped
        manta_small_variants
        manta_small_variants_index
        manta_small_variants_bgzipped
        manta_candidate_variants
        manta_candidate_variants_index
        manta_candidate_variants_bgzipped
        manta_diploid_variants
        manta_diploid_variants_index
        manta_diploid_variants_bgzipped
        gridss_variants
        gridss_variants_index
        gridss_variants_bgzipped
        dysgu_variants
        dysgu_variants_index
        dysgu_variants_bgzipped
        tiddit_variants
        tiddit_ploidy
        tiddit_variants_index
        tiddit_variants_bgzipped
        cutesv_variants
        cutesv_variants_index
        cutesv_variants_bgzipped
        delly_flag
        manta_flag
        gridss_flag
        dysgu_flag
        tiddit_flag
        svaba_flag
        sniffles_flag
        cutesv_flag
        debug_flag
        ch_max_distance_breakpoints
        ch_min_supporting_callers
        ch_account_for_type
        ch_account_for_sv_strands
        ch_estimate_distanced_by_sv_size
        ch_min_sv_size
        ch_taxonomy_id
        ch_outdir

    main:
        view("🚀 Starting REPORT_GENERATION workflow"
        )
        if (debug_flag) {
            genome_file.view            { "DEBUG -> genome_file: ${it.baseName}" }
            delly_variants.view         { "DEBUG -> delly_variants: ${it}" }
            delly_variants_index.view   { "DEBUG -> delly_variants_index: ${it}" }
            gridss_variants.view        { "DEBUG -> gridss_variants: ${it}" }
            tiddit_variants.view        { "DEBUG -> tiddit_variants: ${it}" }
            tiddit_ploidy.view          { "DEBUG -> tiddit_ploidy: ${it}" }
            cutesv_variants.view        { "DEBUG -> cutesv_variants: ${it}" }
        }

        // STEP 2: Prepare VCF files for SURVIVOR_MERGE
        // Collect all VCF files to be merged and materialize them into a list of file paths

        all_vcf_paths = channel.value([])
        all_tbi_paths = channel.value([])
        all_vcf_bgzipped_paths = channel.value([])
        all_meta = channel.value([])

        if (delly_flag) {
            path_delly_variants             = delly_variants.map { it[1] }
            path_delly_variants_bgzipped    = delly_variants_bgzipped.map { it[1] }
            path_delly_variants_index       = delly_variants_index.map { it[1] }
            meta_delly                      = delly_variants.map { it[0] }

            all_vcf_paths = all_vcf_paths.concat(path_delly_variants)
            all_vcf_bgzipped_paths = all_vcf_bgzipped_paths.concat(path_delly_variants_bgzipped)
            all_tbi_paths = all_tbi_paths.concat(path_delly_variants_index)
            all_meta = all_meta.concat(meta_delly)
        }

        if (manta_flag) {
            path_manta_small_variants             = manta_small_variants.map { it[1] }
            path_manta_small_variants_bgzipped    = manta_small_variants_bgzipped.map { it[1] }
            path_manta_small_variants_index       = manta_small_variants_index.map { it[1] }
            meta_manta_small                      = manta_small_variants.map { it[0] }

            all_vcf_paths = all_vcf_paths.concat(path_manta_small_variants)
            all_vcf_bgzipped_paths = all_vcf_bgzipped_paths.concat(path_manta_small_variants_bgzipped)
            all_tbi_paths = all_tbi_paths.concat(path_manta_small_variants_index)
            all_meta = all_meta.concat(meta_manta_small)

            path_manta_candidate_variants             = manta_candidate_variants.map { it[1] }
            path_manta_candidate_variants_bgzipped    = manta_candidate_variants_bgzipped.map { it[1] }
            path_manta_candidate_variants_index       = manta_candidate_variants_index.map { it[1] }
            meta_manta_candidate                      = manta_candidate_variants.map { it[0] }

            all_vcf_paths = all_vcf_paths.concat(path_manta_candidate_variants)
            all_vcf_bgzipped_paths = all_vcf_bgzipped_paths.concat(path_manta_candidate_variants_bgzipped)
            all_tbi_paths = all_tbi_paths.concat(path_manta_candidate_variants_index)
            all_meta = all_meta.concat(meta_manta_candidate)

            path_manta_diploid_variants             = manta_diploid_variants.map { it[1] }
            path_manta_diploid_variants_bgzipped    = manta_diploid_variants_bgzipped.map { it[1] }
            path_manta_diploid_variants_index       = manta_diploid_variants_index.map { it[1] }
            meta_manta_diploid                      = manta_diploid_variants.map { it[0] }

            all_vcf_paths = all_vcf_paths.concat(path_manta_diploid_variants)
            all_vcf_bgzipped_paths = all_vcf_bgzipped_paths.concat(path_manta_diploid_variants_bgzipped)
            all_tbi_paths = all_tbi_paths.concat(path_manta_diploid_variants_index)
            all_meta = all_meta.concat(meta_manta_diploid)
        }

        if (gridss_flag) {
            path_gridss_variants            = gridss_variants.map { it[1] }
            path_gridss_variants_bgzipped   = gridss_variants_bgzipped.map { it[1] }
            path_gridss_variants_index      = gridss_variants_index.map { it[1] }
            meta_gridss                     = gridss_variants.map { it[0] }

            all_vcf_paths = all_vcf_paths.concat(path_gridss_variants)
            all_vcf_bgzipped_paths = all_vcf_bgzipped_paths.concat(path_gridss_variants_bgzipped)
            all_tbi_paths = all_tbi_paths.concat(path_gridss_variants_index)
            all_meta = all_meta.concat(meta_gridss)
        }

        if (dysgu_flag) {
            path_dysgu_variants            = dysgu_variants.map { it[1] }
            path_dysgu_variants_bgzipped   = dysgu_variants_bgzipped.map { it[1] }
            path_dysgu_variants_index      = dysgu_variants_index.map { it[1] }
            meta_dysgu                     = dysgu_variants.map { it[0] }

            all_vcf_paths = all_vcf_paths.concat(path_dysgu_variants)
            all_vcf_bgzipped_paths = all_vcf_bgzipped_paths.concat(path_dysgu_variants_bgzipped)
            all_tbi_paths = all_tbi_paths.concat(path_dysgu_variants_index)
            all_meta = all_meta.concat(meta_dysgu)
        }

        if (tiddit_flag) {
            path_tiddit_variants            = tiddit_variants.map { it[1] }
            path_tiddit_variants_bgzipped   = tiddit_variants_bgzipped.map { it[1] }
            path_tiddit_variants_index      = tiddit_variants_index.map { it[1] }
            meta_tiddit                     = tiddit_variants.map { it[0] }

            all_vcf_paths = all_vcf_paths.concat(path_tiddit_variants)
            all_vcf_bgzipped_paths = all_vcf_bgzipped_paths.concat(path_tiddit_variants_bgzipped)
            all_tbi_paths = all_tbi_paths.concat(path_tiddit_variants_index)
            all_meta = all_meta.concat(meta_tiddit)
        }

        if (cutesv_flag) {
            path_cutesv_variants            = cutesv_variants.map { it[1] }
            path_cutesv_variants_bgzipped   = cutesv_variants_bgzipped.map { it[1] }
            path_cutesv_variants_index      = cutesv_variants_index.map { it[1] }
            meta_cutesv                     = cutesv_variants.map { it[0] }

            all_vcf_paths = all_vcf_paths.concat(path_cutesv_variants)
            all_vcf_bgzipped_paths = all_vcf_bgzipped_paths.concat(path_cutesv_variants_bgzipped)
            all_tbi_paths = all_tbi_paths.concat(path_cutesv_variants_index)
            all_meta = all_meta.concat(meta_cutesv)
        }
\
        all_vcf_paths = all_vcf_paths.filter { it }.toList()
        all_vcf_bgzipped_paths = all_vcf_bgzipped_paths.filter { it }.toList()
        all_tbi_paths = all_tbi_paths.filter { it }.toList()
        all_meta = all_meta.filter { it }.toList()

        meta_custom = Channel.of([id: "survivor_merge"])

        survivor_input = meta_custom.combine(all_vcf_paths.toList())


            // meta_custom.view    { "DEBUG -> meta_custom: ${it}" }
            // all_vcf_paths.view  { "DEBUG -> all_vcf_paths: ${it}" }
            // survivor_input.view { "DEBUG -> survivor_input: ${it}" }

        merged_variants = SURVIVOR_MERGE(
            survivor_input,
            ch_max_distance_breakpoints,
            ch_min_supporting_callers,
            ch_account_for_type,
            ch_account_for_sv_strands,
            ch_estimate_distanced_by_sv_size,
            ch_min_sv_size
        ).vcf

        stats_out = SURVIVOR_STATS(
            merged_variants,
            -1,
            -1,
            -1
        )

        if (debug_flag) {
            stats_out.stats.view { "DEBUG -> stats_out: ${it}" }
        }

        // STEP 3: Merge VCF files with BCFTOOLS
        bfcmerge_meta = Channel.of([id: "bfcmerge_merge"])

        bcfmerge_input = bfcmerge_meta.combine(all_vcf_bgzipped_paths.toList()).combine(all_tbi_paths.toList())

        if (debug_flag) {
            bcfmerge_input.view { "DEBUG -> bcfmerge_input: ${it}" }
        }

        bcftool_report = BCFTOOLS_MERGE(
            bcfmerge_input,
            [[], []],
            [[], []],
            [[], []]
        )

        // STEP 4: Generate HTML report
        html_report = VCF_REPGEN(
            bcftool_report.vcf.map { it[1] },
            all_vcf_paths.flatten(),
            merged_variants.map { it[1] },
            stats_out.stats.map { it[1] },
            ch_taxonomy_id,
            ch_outdir
        )

        // nf-core modules install annotsv/annotsv
        // nf-core modules install svanalyzer/svbenchmark


    emit:
        html_index = html_report.html_index
        html_merged = html_report.html_merged
        html_survivor = html_report.html_survivor

}
