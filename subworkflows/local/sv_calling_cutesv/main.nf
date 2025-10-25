/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: SV_CALLING_CUTESV
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Fast and accurate long-read structural variant detection using cuteSV
    Read type: long-read
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Description:
        Detects structural variants from long-read sequencing data (PacBio/Nanopore)
        using cuteSV's efficient signature-clustering algorithm. Provides fast and accurate
        SV calling optimized for long-read data with low computational requirements.

    Processing Steps:
        1. CUTESV - Calls SVs from long-read BAM files
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

include { CUTESV            } from '../../../modules/nf-core/cutesv/main'
include { SVYNC             } from '../../../modules/nf-core/svync/main'
include { SAMPLE_REHEADER   } from '../../../modules/local/sample_regen/main.nf'

workflow SV_CALLING_CUTESV {
    take:
        bam_inputs
        reference_genome_bgzipped

    main:
        ch_versions = channel.empty()
        name_cutesv = "cutesv"

        CUTESV(
            bam_inputs.map { meta, bam, bai ->
                tuple(meta + [id: "${meta.id}-${name_cutesv}"], bam, bai)
            },
            reference_genome_bgzipped
        )

        SVYNC(
            CUTESV.out.vcf
                .map { meta, vcf ->
                    tuple(meta + [id: "${meta.id}-svync"], vcf, [])
                }
                .combine(
                    channel.value(file("${projectDir}/assets/svync/${name_cutesv}.yaml"))
                )
        )

        SAMPLE_REHEADER(
            SVYNC.out.vcf,
            name_cutesv
        )

        ch_versions = ch_versions.mix(CUTESV.out.versions)
        ch_versions = ch_versions.mix(SVYNC.out.versions)
        ch_versions = ch_versions.mix(SAMPLE_REHEADER.out.versions)

    emit:
        vcf = SAMPLE_REHEADER.out.vcf
        vcfgz =  SAMPLE_REHEADER.out.vcfgz
        tbi = SAMPLE_REHEADER.out.tbi
        csi = SAMPLE_REHEADER.out.csi
        versions = ch_versions
}
