/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: SV_CALLING_DELLY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Structural variant detection using Delly for paired-end sequencing data
    Read type: short-read
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Description:
        Detects structural variants including deletions, duplications, inversions, and
        translocations using Delly's split-read and paired-end mapping approach. The
        workflow synchronizes and standardizes VCF output for downstream analysis.

    Processing Steps:
        1. DELLY_CALL - Detects SVs from aligned BAM files
        2. SVYNC - Synchronizes and refines SV calls to standard format
        3. SAMPLE_REHEADER - Reheaders and renames the output VCF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Input Channels (take):
        bam_bai                       [meta, bam, bai]      Aligned BAM files with index
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

include { DELLY_CALL        } from '../../../modules/nf-core/delly/call/main.nf'
include { SVYNC             } from '../../../modules/nf-core/svync/main'
include { SAMPLE_REHEADER   } from '../../../modules/local/sample_regen/main.nf'

workflow SV_CALLING_DELLY {

    take:
        bam_bai
        reference_genome_unzipped
        reference_genome_faidx

    main:
        ch_versions = Channel.empty()
        name_delly = "delly"

        DELLY_CALL(
            bam_bai.map { meta, bam, bai ->
                tuple(meta + [id: "${meta.id}-${name_delly}"], bam, bai, [], [], [])
            },
            reference_genome_unzipped,
            reference_genome_faidx
        )

        SVYNC(
            DELLY_CALL.out.bcf
                .join(DELLY_CALL.out.csi, by: 0)
                .map { meta, vcf, tbi ->
                    tuple(meta + [id: "${meta.id}-svync"], vcf, tbi)
                }
                .combine(
                    Channel.value(file("${projectDir}/assets/svync/${name_delly}.yaml"))
                )
        )

        SAMPLE_REHEADER(
            SVYNC.out.vcf,
            name_delly
        )

        ch_versions = ch_versions.mix(DELLY_CALL.out.versions.first())
        ch_versions = ch_versions.mix(SVYNC.out.versions.first())
        ch_versions = ch_versions.mix(SAMPLE_REHEADER.out.versions.first())

    emit:
        vcf = SAMPLE_REHEADER.out.vcf
        vcfgz =  SAMPLE_REHEADER.out.vcfgz
        tbi = SAMPLE_REHEADER.out.tbi
        csi = SAMPLE_REHEADER.out.csi
        versions = ch_versions
}
