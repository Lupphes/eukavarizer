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
        bcfmerge_vcf              tuple(meta, vcf.gz)    BCFtools merged and filtered VCF
        bcfmerge_tbi              tuple(meta, tbi)       BCFtools VCF index
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

workflow SV_UNIFICATION {

    take:
        vcf_list
        vcfgz_list
        tbi_list
        reference_genome_bgzipped

    main:
        ch_versions = Channel.empty()

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

        ch_versions = ch_versions.mix(SURVIVOR_MERGE.out.versions.first())
        ch_versions = ch_versions.mix(SURVIVOR_FILTER.out.versions.first())
        ch_versions = ch_versions.mix(SURVIVOR_STATS.out.versions.first())
        ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions.first())
        ch_versions = ch_versions.mix(BCFTOOLS_FILTER.out.versions.first())
        ch_versions = ch_versions.mix(BCFTOOLS_SORT.out.versions.first())
        ch_versions = ch_versions.mix(BCFTOOLS_STATS.out.versions.first())

    emit:
        survivor_vcf    = SURVIVOR_FILTER.out.vcf
        survivor_stats  = SURVIVOR_STATS.out.stats
        bcfmerge_vcf    = BCFTOOLS_SORT.out.vcf
        bcfmerge_tbi    = BCFTOOLS_SORT.out.tbi
        bcfmerge_stats  = BCFTOOLS_STATS.out.stats
        versions        = ch_versions
}
