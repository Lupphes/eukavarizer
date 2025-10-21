/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: SV_CALLING_DYSGU
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Machine learning-based structural variant detection using Dysgu
    Read type: short-read
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Description:
        Detects structural variants using Dysgu's machine learning model trained on
        diverse datasets. Provides accurate SV calling for short-read sequencing data
        through advanced algorithmic approaches and ML-based classification.

    Processing Steps:
        1. DYSGU - Detects SVs using machine learning from BAM files
        2. SVYNC - Synchronizes and refines SV calls to standard format
        3. SAMPLE_REHEADER - Reheaders and renames the output VCF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Input Channels (take):
        bam_inputs                    [meta, bam, bai]      Aligned BAM files with index
        reference_genome_unzipped     [meta, fasta]         Reference genome FASTA
        reference_genome_faidx        [meta, fai]           Reference genome index
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

include { DYSGU             } from '../../../modules/nf-core/dysgu/main'
include { SVYNC             } from '../../../modules/nf-core/svync/main'
include { SAMPLE_REHEADER   } from '../../../modules/local/sample_regen/main.nf'

workflow SV_CALLING_DYSGU {
    take:
        bam_inputs
        reference_genome_unzipped
        reference_genome_faidx

    main:
        ch_versions = Channel.empty()
        name_dysgu = "dysgu"

        DYSGU(
            bam_inputs
            .map { meta, bam, bai ->
                tuple(meta + [id: "${meta.id}-${name_dysgu}"], bam, bai)
            },
            reference_genome_unzipped
            .combine(reference_genome_faidx)
            .map { meta_fasta, fasta, _meta_fai, fai -> [meta_fasta, fasta, fai] }
            .collect(),
        )

        SVYNC(
            DYSGU.out.vcf
                .join(DYSGU.out.tbi, by: 0)
                .map { meta, vcf, tbi ->
                    tuple(meta + [id: "${meta.id}-svync"], vcf, tbi)
                }
                .combine(
                    Channel.value(file("${projectDir}/assets/svync/${name_dysgu}.yaml"))
                )
        )

        SAMPLE_REHEADER(
            SVYNC.out.vcf,
            name_dysgu
        )

        ch_versions = ch_versions.mix(DYSGU.out.versions.first())
        ch_versions = ch_versions.mix(SVYNC.out.versions.first())
        ch_versions = ch_versions.mix(SAMPLE_REHEADER.out.versions.first())

    emit:
        vcf = SAMPLE_REHEADER.out.vcf
        vcfgz =  SAMPLE_REHEADER.out.vcfgz
        tbi = SAMPLE_REHEADER.out.tbi
        csi = SAMPLE_REHEADER.out.csi
        versions = ch_versions

}
