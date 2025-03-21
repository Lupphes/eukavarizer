/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { VARIFY            } from '../../../modules/local/varify/main'
// include { DUPHOLD               } from '../../modules/nf-core/duphold/main.nf'

workflow REPORT_GENERATION {

    take:
        ch_taxonomy_id
        ch_outdir
        vcf_list
        survivor_vcf
        survivor_stats
        bcfmerge_vcf

    main:
        VARIFY(
            bcfmerge_vcf.map { it[1] },
            vcf_list.filter { it }
            .map { it[1] }
            .toList().flatten(),
            survivor_vcf.map { it[1] },
            survivor_stats.map { it[1] },
            ch_taxonomy_id,
            ch_outdir
        )

    emit:
        html_index      = VARIFY.out.html_index
        html_merged     = VARIFY.out.html_merged
        html_survivor   = VARIFY.out.html_survivor




        // first = fastq_bam
        //     .join(fastq_bam_indexes, by: 0)
        //     .map { meta, bam, bai -> tuple(meta, bam, bai, vcf, [], []) }

        // DUPHOLD(
        //     first,
        //     fastq_bam_indexes.map { _meta, ref -> ref },
        //     reference_genome_bgzipped.map { _meta, ref -> ref }
        // )
}
