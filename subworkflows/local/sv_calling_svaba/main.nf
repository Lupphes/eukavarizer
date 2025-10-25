/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: SV_CALLING_SVABA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Local assembly-based SV and indel detection using SvABA for Illumina data
    Read type: short-read
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Description:
        Detects structural variants and indels using SvABA's local assembly approach,
        designed for Illumina paired-end data. Provides sensitive detection of complex
        rearrangements through local micro-assembly and realignment.

    Processing Steps:
        1. SVABA - Calls SVs and indels using local assembly from BAM files
        2. SVABA_ANNOTATE (optional) - Annotates SVs with additional information
        3. SVYNC - Synchronizes and refines SV calls to standard format
        4. SAMPLE_REHEADER - Reheaders and renames the output VCF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Input Channels (take):
        bam_inputs                    [meta, bam, bai]      Aligned BAM files with index
        reference_genome_unzipped     [meta, fasta]         Reference genome FASTA
        reference_genome_faidx        [meta, fai]           Reference genome index
        reference_genome_bwa_index    [meta, index_dir]     BWA reference index
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

include { SVABA             } from '../../../modules/nf-core/svaba/main'
include { SVABA_ANNOTATE    } from '../../../modules/local/svaba/annotate/main.nf'
include { SVYNC             } from '../../../modules/nf-core/svync/main'
include { SAMPLE_REHEADER   } from '../../../modules/local/sample_regen/main.nf'

workflow SV_CALLING_SVABA {
    take:
        bam_inputs
        reference_genome_unzipped
        reference_genome_faidx
        reference_genome_bwa_index

    main:
        ch_versions = channel.empty()
        name_svaba = "svaba"

        SVABA(
            // SVABA requires Illumina paired-end reads, does not support minimap2 & bwamem2
            bam_inputs.filter { meta, _bam, _bai ->
                !(params.bwamem2) &&
                (
                    (meta.platform == 'illumina') ||
                    (!meta.platform && meta.median_bp <= params.long_read_threshold)
                )
            }.map { meta, bam, bai ->
                tuple(meta + [id: "${meta.id}-${name_svaba}"], bam, bai, [], [])
            },
            reference_genome_unzipped,
            reference_genome_faidx,
            reference_genome_bwa_index,
            [[],[]],
            [[],[]],
            [[],[]],
        )

        if (params.svaba_annotate) {
            svaba_annotate = SVABA_ANNOTATE(
                SVABA.out.sv
            ).vcf
        } else {
            svaba_annotate = SVABA.out.sv
        }

        SVYNC(
            svaba_annotate
                .map { meta, vcf ->
                    tuple(meta + [id: "${meta.id}-${name_svaba}-svync"], vcf, [])
                }
                .combine(
                    channel.value(file("${projectDir}/assets/svync/${name_svaba}.yaml"))
                )
        )

        SAMPLE_REHEADER(
            SVYNC.out.vcf,
            name_svaba
        )

        ch_versions = ch_versions.mix(SVABA.out.versions)
        if (params.svaba_annotate) {
            ch_versions = ch_versions.mix(SVABA_ANNOTATE.out.versions)
        }
        ch_versions = ch_versions.mix(SVYNC.out.versions)
        ch_versions = ch_versions.mix(SAMPLE_REHEADER.out.versions)

    emit:
        vcf = SAMPLE_REHEADER.out.vcf
        vcfgz =  SAMPLE_REHEADER.out.vcfgz
        tbi = SAMPLE_REHEADER.out.tbi
        csi = SAMPLE_REHEADER.out.csi
        versions = ch_versions
}
