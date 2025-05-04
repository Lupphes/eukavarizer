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
    - `survivor_vcf`     – Filtered and merged VCF from SURVIVOR.
    - `survivor_stats`   – Summary statistics from SURVIVOR.
    - `bcfmerge_vcf`     – Filtered merged VCF from BCFTOOLS.
    - `bcfmerge_tbi`     – Tabix index for the BCFTOOLS merged VCF.
    - `bcfmerge_stats`   – Summary statistics from BCFTOOLS.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SURVIVOR_MERGE        } from '../../../modules/nf-core/survivor/merge/main'
include { SURVIVOR_STATS        } from '../../../modules/nf-core/survivor/stats/main'
include { SURVIVOR_FILTER       } from '../../../modules/nf-core/survivor/filter/main'

include { BCFTOOLS_CONCAT       } from '../../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_SORT         } from '../../../modules/nf-core/bcftools/sort/main'
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
        vcf_list_cleaned = vcf_list
            .filter { it != null }
            .map { it[1] }

        vcfgz_list_cleaned = vcfgz_list
            .filter { it != null }
            .map { it[1] }

        tbi_list_cleaned = tbi_list
            .filter { it != null }
            .map { it[1] }


        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            SURVIVOR MERGE
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        survivor_input = vcf_list_cleaned
            .collect().unique()
            .map { vcf_files ->
                tuple([id: "survivor_merge"], vcf_files)
            }

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

        bcfmerge_input = vcfgz_list_cleaned.collect().unique().collect(flat: false)
            .combine(tbi_list_cleaned.collect().unique().collect(flat: false))
            .filter { vcfgz_files, tbi_files -> vcfgz_files && tbi_files }
            .map { vcfgz_files, tbi_files ->
                tuple([id: "bcftools_concat"], vcfgz_files, tbi_files)
            }

        BCFTOOLS_CONCAT(
            bcfmerge_input
        )

        BCFTOOLS_FILTER(
            BCFTOOLS_CONCAT.out.vcf
            .join(BCFTOOLS_CONCAT.out.tbi, by: 0)
        )

        BCFTOOLS_SORT(
            BCFTOOLS_FILTER.out.vcf
        )

        BCFTOOLS_STATS(
            BCFTOOLS_SORT.out.vcf
            .join(BCFTOOLS_SORT.out.tbi, by: 0),
            [[id: "regions"], []],
            [[id: "targets"], []],
            [[id: "samples"], []],
            [[id: "exons"], []],
            reference_genome_bgzipped
        )

    emit:
        survivor_vcf    = SURVIVOR_FILTER.out.vcf
        survivor_stats  = SURVIVOR_STATS.out.stats
        bcfmerge_vcf    = BCFTOOLS_SORT.out.vcf
        bcfmerge_tbi    = BCFTOOLS_SORT.out.tbi
        bcfmerge_stats  = BCFTOOLS_STATS.out.stats
}
