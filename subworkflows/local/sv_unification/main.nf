/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: SV_UNIFICATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Merges and unifies structural variant calls from multiple callers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Description:
        This subworkflow combines structural variant calls from different tools using
        two complementary approaches: SURVIVOR for consensus-based merging and BCFtools
        for sample-based merging. Both outputs are filtered and annotated with statistics.

    Processing Steps:
        1. Merge VCF files using SURVIVOR based on breakpoint distance and caller support
        2. Filter SURVIVOR merged VCF by size and quality thresholds
        3. Generate statistics for SURVIVOR merged variants
        4. Merge VCF files using BCFtools for sample-level consolidation
        5. Filter BCFtools merged VCF and generate statistics
        6. (Optional) Collapse redundant variants using Truvari
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Input Channels (take):
        vcf_list                  tuple(meta, vcf)       Individual caller VCF files
        vcfgz_list                tuple(meta, vcf.gz)    Compressed caller VCF files
        tbi_list                  tuple(meta, tbi)       VCF index files
        reference_genome_bgzipped tuple(meta, fasta.gz)  Compressed reference genome
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Output Channels (emit):
        survivor_vcf              tuple(meta, vcf)       SURVIVOR merged and filtered VCF
        survivor_stats            tuple(meta, stats)     SURVIVOR statistics
        concat_vcf                tuple(meta, vcf.gz)    BCFtools concat with truvari
        concat_tbi                tuple(meta, tbi)       BCFtools VCF index
        bcfmerge_stats            tuple(meta, stats)     BCFtools statistics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Author:   Ondřej Sloup (Lupphes)
    Contact:  ondrej.sloup@protonmail.com
    GitHub:   @Lupphes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SURVIVOR_MERGE        } from '../../../modules/nf-core/survivor/merge/main'
include { SURVIVOR_STATS        } from '../../../modules/nf-core/survivor/stats/main'
include { SURVIVOR_FILTER       } from '../../../modules/nf-core/survivor/filter/main'

include { BCFTOOLS_CONCAT       } from '../../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_SORT         } from '../../../modules/nf-core/bcftools/sort/main'
include { BCFTOOLS_STATS        } from '../../../modules/nf-core/bcftools/stats/main'
include { BCFTOOLS_FILTER       } from '../../../modules/nf-core/bcftools/filter/main'

include { TRUVARI_COLLAPSE      } from '../../../modules/local/truvari/collapse/main'
include { BCFTOOLS_SORT as BCFTOOLS_SORT_TRUVARI    } from '../../../modules/nf-core/bcftools/sort/main'
include { TABIX_TABIX           } from '../../../modules/nf-core/tabix/tabix/main'

workflow SV_UNIFICATION {

    take:
        vcf_list
        vcfgz_list
        tbi_list
        reference_genome_bgzipped
        reference_genome_faidx

    main:
        ch_versions = channel.empty()

        vcf_list_cleaned = vcf_list
            .filter { tuple -> tuple != null }
            .map { tuple -> tuple[1] }

        vcfgz_list_cleaned = vcfgz_list
            .filter { tuple -> tuple != null }
            .map { tuple -> tuple[1] }

        tbi_list_cleaned = tbi_list
            .filter { tuple -> tuple != null }
            .map { tuple -> tuple[1] }


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
        
        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            SURVIVOR STATS
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

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

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            TRUVARI COLLAPSE
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        
        if (params.truvari_collapse) {
            TRUVARI_COLLAPSE(
                BCFTOOLS_SORT.out.vcf
                .join(BCFTOOLS_SORT.out.tbi, by: 0),
                reference_genome_bgzipped,
                reference_genome_faidx
            )

            BCFTOOLS_SORT_TRUVARI(
                TRUVARI_COLLAPSE.out.collapsed_vcf
            )

            ch_versions = ch_versions.mix(TRUVARI_COLLAPSE.out.versions)
            ch_versions = ch_versions.mix(BCFTOOLS_SORT_TRUVARI.out.versions)
            concat_vcf = BCFTOOLS_SORT_TRUVARI.out.vcf
            concat_tbi = BCFTOOLS_SORT_TRUVARI.out.tbi
        } else {
            concat_vcf = BCFTOOLS_SORT.out.vcf
            concat_tbi = BCFTOOLS_SORT.out.tbi
        }

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            BCFTOOLS STATS
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        BCFTOOLS_STATS(
            concat_vcf
            .join(concat_tbi, by: 0),
            [[id: "regions"], []],
            [[id: "targets"], []],
            [[id: "samples"], []],
            [[id: "exons"], []],
            reference_genome_bgzipped
        )

        ch_versions = ch_versions.mix(SURVIVOR_MERGE.out.versions)
        ch_versions = ch_versions.mix(SURVIVOR_FILTER.out.versions)
        ch_versions = ch_versions.mix(SURVIVOR_STATS.out.versions)
        ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions)
        ch_versions = ch_versions.mix(BCFTOOLS_FILTER.out.versions)
        ch_versions = ch_versions.mix(BCFTOOLS_SORT.out.versions)
        ch_versions = ch_versions.mix(BCFTOOLS_STATS.out.versions)

    emit:
        survivor_vcf    = SURVIVOR_FILTER.out.vcf
        survivor_stats  = SURVIVOR_STATS.out.stats
        concat_vcf      = concat_vcf
        concat_tbi      = concat_tbi
        bcfmerge_stats  = BCFTOOLS_STATS.out.stats
        versions        = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
