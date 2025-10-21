/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: SV_CALLING_GRIDSS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    High-precision structural variant detection using local assembly with GRIDSS
    Read type: short-read
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Description:
        Detects structural variants using GRIDSS's local assembly approach for germline
        and somatic samples. Provides highly accurate breakpoint detection through local
        assembly and split-read analysis, with optional variant annotation.

    Processing Steps:
        1. GRIDSS_GRIDSS - Detects SVs using local assembly from BAM files
        2. GRIDSS_ANNOTATE (optional) - Annotates SVs with additional information
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

include { GRIDSS_GRIDSS     } from '../../../modules/nf-core/gridss/gridss/main'
include { GRIDSS_ANNOTATE   } from '../../../modules/local/gridss/annotate/main.nf'
include { SVYNC             } from '../../../modules/nf-core/svync/main'
include { SAMPLE_REHEADER   } from '../../../modules/local/sample_regen/main.nf'

workflow SV_CALLING_GRIDSS {
    take:
        bam_inputs
        reference_genome_unzipped
        reference_genome_faidx
        reference_genome_bwa_index

    main:
        ch_versions = Channel.empty()
        name_gridss = "gridss"

        GRIDSS_GRIDSS(
            bam_inputs
                // Don't use on long reads (PacBio/Nanopore raw reads). GRIDSS does not support MINIMAP2.
                .filter { meta, _bam, _bai ->
                    !params.minimap2_flag && (
                        (meta.platform == 'illumina') ||
                        (!meta.platform && meta.median_bp <= params.long_read_threshold)
                    )
                }
                .map { meta, bam, _bai -> tuple(meta + [id: "${meta.id}-${name_gridss}"], bam) },
            reference_genome_unzipped,
            reference_genome_faidx,
            // reference_genome_minimap_index -> GRIDSS does not support MINIMAP2
            reference_genome_bwa_index
        )

        if (params.gridss_annotate) {
            gridss_annotate = GRIDSS_ANNOTATE(
                GRIDSS_GRIDSS.out.vcf
            ).vcf
        } else {
            gridss_annotate = GRIDSS_GRIDSS.out.vcf
        }

        SVYNC(
            gridss_annotate
                .map { meta, vcf ->
                    tuple(meta + [id: "${meta.id}-svync"], vcf, [])
                }
                .combine(
                    Channel.value(file("${projectDir}/assets/svync/${name_gridss}.yaml"))
                )
        )

        SAMPLE_REHEADER(
            SVYNC.out.vcf,
            name_gridss
        )

        ch_versions = ch_versions.mix(GRIDSS_GRIDSS.out.versions.first())
        if (params.gridss_annotate) {
            ch_versions = ch_versions.mix(GRIDSS_ANNOTATE.out.versions.first())
        }
        ch_versions = ch_versions.mix(SVYNC.out.versions.first())
        ch_versions = ch_versions.mix(SAMPLE_REHEADER.out.versions.first())

    emit:
        vcf = SAMPLE_REHEADER.out.vcf
        vcfgz =  SAMPLE_REHEADER.out.vcfgz
        tbi = SAMPLE_REHEADER.out.tbi
        csi = SAMPLE_REHEADER.out.csi
        versions = ch_versions
}
