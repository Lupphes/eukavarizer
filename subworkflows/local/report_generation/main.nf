// include { DUPHOLD               } from '../../modules/nf-core/duphold/main.nf'

workflow REPORT_GENERATION {


}


        // first = fastq_bam
        //     .join(fastq_bam_indexes, by: 0)
        //     .map { meta, bam, bai -> tuple(meta, bam, bai, vcf, [], []) }

        // DUPHOLD(
        //     first,
        //     fastq_bam_indexes.map { _meta, ref -> ref },
        //     reference_genome_bgzipped.map { _meta, ref -> ref }
        // )
