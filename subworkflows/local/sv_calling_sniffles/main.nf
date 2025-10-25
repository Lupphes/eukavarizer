/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: SV_CALLING_SNIFFLES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Long-read structural variant detection using Sniffles for PacBio and Nanopore data
    Read type: long-read
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Description:
        Detects structural variants from long-read sequencing data (PacBio/Nanopore)
        using Sniffles. Leverages split-read signatures and coverage analysis to identify
        SVs with high sensitivity and accuracy from long-read alignments.

    Processing Steps:
        1. SNIFFLES - Calls SVs from long-read BAM files
        2. SVYNC - Synchronizes and refines SV calls to standard format
        3. SAMPLE_REHEADER - Reheaders and renames the output VCF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Input Channels (take):
        bam_inputs                    [meta, bam, bai]      Aligned BAM files with index
        reference_genome_bgzipped     [meta, fasta.gz]      Bgzipped reference genome
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Output Channels (emit):
        vcf                           [meta, vcf]           Reheaded VCF file
        vcfgz                         [meta, vcf.gz]        Gzipped VCF file
        tbi                           [meta, tbi]           Tabix index
        csi                           [meta, csi]           CSI index
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Author:   Ondřej Sloup (Lupphes)
    Contact:  ondrej.sloup@protonmail.com
    GitHub:   @Lupphes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SNIFFLES          } from '../../../modules/nf-core/sniffles/main'
include { SVYNC             } from '../../../modules/nf-core/svync/main'
include { SAMPLE_REHEADER   } from '../../../modules/local/sample_regen/main.nf'

workflow SV_CALLING_SNIFFLES {
    take:
        bam_inputs
        reference_genome_bgzipped

    main:
        ch_versions = channel.empty()
        name_sniffles = "sniffles"

        SNIFFLES(
            bam_inputs
            .map { meta, bam, bai ->
                tuple(meta + [id: "${meta.id}-${name_sniffles}"], bam, bai)
            },
            reference_genome_bgzipped,
            [[],[]],
            true,
            false
        )

        SVYNC(
            SNIFFLES.out.vcf
                .join(SNIFFLES.out.tbi, by: 0)
                .map { meta, vcf, tbi ->
                    tuple(meta + [id: "${meta.id}-svync"], vcf, tbi)
                }
                .combine(
                    channel.value(file("${projectDir}/assets/svync/${name_sniffles}.yaml"))
                )
        )

        SAMPLE_REHEADER(
            SVYNC.out.vcf,
            name_sniffles
        )

        ch_versions = ch_versions.mix(SNIFFLES.out.versions)
        ch_versions = ch_versions.mix(SVYNC.out.versions)
        ch_versions = ch_versions.mix(SAMPLE_REHEADER.out.versions)

    emit:
        vcf = SAMPLE_REHEADER.out.vcf
        vcfgz =  SAMPLE_REHEADER.out.vcfgz
        tbi = SAMPLE_REHEADER.out.tbi
        csi = SAMPLE_REHEADER.out.csi
        versions = ch_versions
}
