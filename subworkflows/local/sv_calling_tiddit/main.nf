/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: SV_CALLING_TIDDIT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Coverage-based structural variant detection using TIDDIT for paired-end reads
    Read type: short-read
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Description:
        Detects structural variants using TIDDIT's coverage-based approach combined with
        paired-end and split-read analysis. Supports both BWA and Minimap2 aligned data
        for paired-end Illumina sequencing.

    Processing Steps:
        1. TIDDIT_SV (BWA/MAP) - Detects SVs from paired-end BAM files
        2. SVYNC - Synchronizes and refines SV calls to standard format
        3. SAMPLE_REHEADER - Reheaders and renames the output VCF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Input Channels (take):
        bam_inputs                    [meta, bam, bai]      Aligned BAM files with index
        reference_genome_bgzipped     [meta, fasta.gz]      Bgzipped reference genome
        reference_genome_bwa_index    [meta, index_dir]     BWA reference index
        reference_genome_minimap_index [meta, mmi]          Minimap2 reference index
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

include { TIDDIT_SV as TIDDIT_MAP   } from '../../../modules/nf-core/tiddit/sv/main'
include { TIDDIT_SV as TIDDIT_BWA   } from '../../../modules/nf-core/tiddit/sv/main'
include { SVYNC                     } from '../../../modules/nf-core/svync/main'
include { SAMPLE_REHEADER           } from '../../../modules/local/sample_regen/main.nf'

workflow SV_CALLING_TIDDIT {
    take:
        bam_inputs
        reference_genome_bgzipped
        reference_genome_bwa_index
        reference_genome_minimap_index

    main:
        ch_versions = channel.empty()
        name_tiddit = "tiddit"

        TIDDIT_BWA(
            bam_inputs
            // Only keep paired-end reads as TIDDIT does not support unpaired reads or long reads
            .filter { meta, _bam, _bai ->
                !meta.single_end
            }
            // Use BWA for short reads
            .filter { meta, _bam, _bai ->
                !params.minimap2_flag && (
                    (meta.platform && meta.platform == 'illumina') ||
                    (!meta.platform && meta.median_bp <= params.long_read_threshold)
                )
            }.map { meta, bam, bai ->
                tuple(meta + [id: "${meta.id}-${name_tiddit}"], bam, bai)
            },
            reference_genome_bgzipped,
            reference_genome_bwa_index
        )

        TIDDIT_MAP(
            bam_inputs
            // Only keep paired-end reads as TIDDIT does not support unpaired reads or long reads
            .filter { meta, _bam, _bai ->
                !meta.single_end
            }
            .filter { meta, _bam, _bai ->
                params.minimap2_flag && (
                    (meta.platform && meta.platform == 'illumina') ||
                    (!meta.platform && meta.median_bp <= params.long_read_threshold)
                )
            }.map { meta, bam, bai ->
                tuple(meta + [id: "${meta.id}-${name_tiddit}"], bam, bai)
            },
            reference_genome_bgzipped,
            reference_genome_minimap_index
        )

        tiddit_result = TIDDIT_BWA.out.vcf.mix(TIDDIT_MAP.out.vcf)

        SVYNC(
            tiddit_result
                .map { meta, vcf ->
                    tuple(meta + [id: "${meta.id}-svync"], vcf, [])
                }
                .combine(
                    channel.value(file("${projectDir}/assets/svync/${name_tiddit}.yaml"))
                )
        )

        SAMPLE_REHEADER(
            SVYNC.out.vcf,
            name_tiddit
        )

        ch_versions = ch_versions.mix(TIDDIT_BWA.out.versions)
        ch_versions = ch_versions.mix(TIDDIT_MAP.out.versions)
        ch_versions = ch_versions.mix(SVYNC.out.versions)
        ch_versions = ch_versions.mix(SAMPLE_REHEADER.out.versions)

    emit:
        vcf = SAMPLE_REHEADER.out.vcf
        vcfgz =  SAMPLE_REHEADER.out.vcfgz
        tbi = SAMPLE_REHEADER.out.tbi
        csi = SAMPLE_REHEADER.out.csi
        versions = ch_versions
}
