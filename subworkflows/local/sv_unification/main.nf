/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    STRUCTURAL VARIANT (SV) MERGING WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    This workflow merges structural variant (SV) calls from different tools:
    1. **SURVIVOR_MERGE** – Merges VCF files and filters based on breakpoint distance,
        caller support, SV type, and strand consistency.
    2. **SURVIVOR_STATS** – Generates summary statistics from the merged VCF.
    3. **BCFTOOLS_MERGE** – Combines compressed VCF files into a unified output.

    Outputs:
    - `survivor_vcf` – Unified VCF from SURVIVOR_MERGE
    - `survivor_stats` – Summary statistics
    - `bcfmerge_vcf` – Unified VCF from BCFTOOLS_MERGE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SURVIVOR_MERGE        } from '../../../modules/nf-core/survivor/merge/main'
include { SURVIVOR_STATS        } from '../../../modules/nf-core/survivor/stats/main'
include { SURVIVOR_FILTER       } from '../../../modules/nf-core/survivor/filter/main'

include { BCFTOOLS_MERGE        } from '../../../modules/nf-core/bcftools/merge/main'
include { BCFTOOLS_STATS        } from '../../../modules/nf-core/bcftools/stats/main'
include { BCFTOOLS_FILTER       } from '../../../modules/nf-core/bcftools/filter/main'

workflow SV_UNIFICATION {

    take:
        vcf_list
        vcfgz_list
        tbi_list
        reference_genome_bgzipped

    main:
        vcf_list = vcf_list.filter { it }
            .map { it[1] }
            .toList()

        vcfgz_list = vcfgz_list.filter { it }
            .map { it[1] }
            .toList()

        tbi_list = tbi_list.filter { it }
            .map { it[1] }
            .toList()

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            SURVIVOR MERGE
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        survivor_meta = Channel.value([id: "survivor_merge"])
        survivor_input = survivor_meta.combine(vcf_list.toList())

        SURVIVOR_MERGE(
            survivor_input,
            params.sur_max_distance_breakpoints,
            params.sur_min_supporting_callers,
            params.sur_account_for_type,
            params.sur_account_for_sv_strands,
            params.sur_estimate_distanced_by_sv_size,
            params.sur_min_sv_size
        )

        SURVIVOR_FILTER(
            SURVIVOR_MERGE.out.vcf.map { meta, vcf -> tuple(meta, vcf, []) },
            params.sur_min_sv_size_filter,
            params.sur_max_sv_size_filter,
            params.sur_min_allele_freq_filter,
            params.sur_min_num_reads_filter
        )

        SURVIVOR_STATS(
            SURVIVOR_FILTER.out.vcf,
            -1,
            -1,
            -1
        )

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            MERGE VCF FILES
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        bfcmerge_meta = Channel.value([id: "bfcmerge_merge"])
        bcfmerge_input = bfcmerge_meta.combine(vcfgz_list.toList()).combine(tbi_list.toList())

        BCFTOOLS_MERGE(
            bcfmerge_input,
            [[], []],
            [[], []],
            [[], []]
        )

        vcf_index_bfcmerge = BCFTOOLS_MERGE.out.vcf
            .map { meta, vfc ->
                tuple(meta, vfc, [])
            }

        BCFTOOLS_FILTER(
            vcf_index_bfcmerge
        )

        vcf_index_bfcmerge_filtered = BCFTOOLS_FILTER.out.vcf
            .map { meta, vfc ->
                tuple(meta, vfc, [])
            }

        BCFTOOLS_STATS(
            vcf_index_bfcmerge_filtered,
            [[], []],
            [[], []],
            [[], []],
            [[], []],
            reference_genome_bgzipped
        )

    emit:
        survivor_vcf = SURVIVOR_FILTER.out.vcf
        survivor_stats = SURVIVOR_STATS.out.stats
        bcfmerge_vcf = BCFTOOLS_FILTER.out.vcf
        bcfmerge_tbi = BCFTOOLS_FILTER.out.tbi
        bcfmerge_stats = BCFTOOLS_STATS.out.stats
}
