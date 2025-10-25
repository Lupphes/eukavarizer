/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: SV_CALLING_MANTA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Germline structural variant detection using Manta for paired-end Illumina data
    Read type: short-read
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Description:
        Detects large structural variants including deletions, inversions, and duplications
        from germline samples using Manta. Outputs three separate VCF streams for small
        indels, candidate SVs, and diploid SVs, each synchronized and standardized.

    Processing Steps:
        1. MANTA_GERMLINE - Detects SVs from paired-end BAM files
        2. SVYNC (3x) - Synchronizes small indels, candidate SVs, and diploid SVs
        3. SAMPLE_REHEADER (3x) - Reheaders and renames output VCFs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Input Channels (take):
        bam_inputs                    [meta, bam, bai]      Aligned BAM files with index
        reference_genome_unzipped     [meta, fasta]         Reference genome FASTA
        reference_genome_faidx        [meta, fai]           Reference genome index
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Output Channels (emit):
        small_vcf                     [meta, vcf]           Small indels VCF
        small_vcfgz                   [meta, vcf.gz]        Small indels gzipped VCF
        small_tbi                     [meta, tbi]           Small indels tabix index
        small_csi                     [meta, csi]           Small indels CSI index
        candidate_vcf                 [meta, vcf]           Candidate SV VCF
        candidate_vcfgz               [meta, vcf.gz]        Candidate SV gzipped VCF
        candidate_tbi                 [meta, tbi]           Candidate SV tabix index
        candidate_csi                 [meta, csi]           Candidate SV CSI index
        diploid_vcf                   [meta, vcf]           Diploid SV VCF
        diploid_vcfgz                 [meta, vcf.gz]        Diploid SV gzipped VCF
        diploid_tbi                   [meta, tbi]           Diploid SV tabix index
        diploid_csi                   [meta, csi]           Diploid SV CSI index
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Author:   Ondřej Sloup (Lupphes)
    Contact:  ondrej.sloup@protonmail.com
    GitHub:   @Lupphes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MANTA_GERMLINE                                } from '../../../modules/nf-core/manta/germline/main'

include { SAMPLE_REHEADER as INDELS_SAMPLE_REHEADER     } from '../../../modules/local/sample_regen/main.nf'
include { SAMPLE_REHEADER as SV_SAMPLE_REHEADER         } from '../../../modules/local/sample_regen/main.nf'
include { SAMPLE_REHEADER as DIPLOID_SAMPLE_REHEADER    } from '../../../modules/local/sample_regen/main.nf'

include { SVYNC as INDEL_SVYNC                          } from '../../../modules/nf-core/svync/main'
include { SVYNC as SV_SVYNC                             } from '../../../modules/nf-core/svync/main'
include { SVYNC as DIPLOID_SVYNC                        } from '../../../modules/nf-core/svync/main'


workflow SV_CALLING_MANTA {
    take:
        bam_inputs
        reference_genome_unzipped
        reference_genome_faidx

    main:
        ch_versions = channel.empty()
        name_manta = "manta"
        name_manta_small = "small"
        name_manta_candidate = "candidate"
        name_manta_diploid = "diploid"

        MANTA_GERMLINE(
            bam_inputs
            // Only keep paired-end reads as MANTA does not support unpaired reads or long reads
            .filter { meta, _bam, _bai ->
                !meta.single_end
            }
            .filter { meta, _bam, _bai ->
                (
                    (meta.platform && meta.platform == 'illumina') ||
                    (!meta.platform && meta.median_bp <= params.long_read_threshold)
                )
            }.map { meta, bam, bai ->
                tuple(meta + [id: "${meta.id}-${name_manta}"], bam, bai, [], [])
            },
            reference_genome_unzipped,
            reference_genome_faidx,
            []
        )

        INDEL_SVYNC(
            MANTA_GERMLINE.out.candidate_small_indels_vcf
                .join(MANTA_GERMLINE.out.candidate_small_indels_vcf_tbi, by: 0)
                .map { meta, vcf, tbi ->
                    tuple(meta + [id: "${meta.id}-${name_manta_small}-svync"], vcf, tbi)
                }
                .combine(
                    channel.value(file("${projectDir}/assets/svync/${name_manta}.yaml"))
                )
        )

        SV_SVYNC(
            MANTA_GERMLINE.out.candidate_sv_vcf
                .join(MANTA_GERMLINE.out.candidate_sv_vcf_tbi, by: 0)
                .map { meta, vcf, tbi ->
                    tuple(meta + [id: "${meta.id}-${name_manta_candidate}-svync"], vcf, tbi)
                }
                .combine(
                    channel.value(file("${projectDir}/assets/svync/${name_manta}.yaml"))
                )
        )

        DIPLOID_SVYNC(
            MANTA_GERMLINE.out.diploid_sv_vcf
                .join(MANTA_GERMLINE.out.diploid_sv_vcf_tbi, by: 0)
                .map { meta, vcf, tbi ->
                    tuple(meta + [id: "${meta.id}-${name_manta_diploid}-svync"], vcf, tbi)
                }
                .combine(
                    channel.value(file("${projectDir}/assets/svync/${name_manta}.yaml"))
                )
        )

        INDELS_SAMPLE_REHEADER(
            INDEL_SVYNC.out.vcf,
            name_manta
        )

        SV_SAMPLE_REHEADER(
            SV_SVYNC.out.vcf,
            name_manta
        )

        DIPLOID_SAMPLE_REHEADER(
            DIPLOID_SVYNC.out.vcf,
            name_manta
        )

        ch_versions = ch_versions.mix(MANTA_GERMLINE.out.versions)
        ch_versions = ch_versions.mix(INDEL_SVYNC.out.versions)
        ch_versions = ch_versions.mix(SV_SVYNC.out.versions)
        ch_versions = ch_versions.mix(DIPLOID_SVYNC.out.versions)
        ch_versions = ch_versions.mix(INDELS_SAMPLE_REHEADER.out.versions)
        ch_versions = ch_versions.mix(SV_SAMPLE_REHEADER.out.versions)
        ch_versions = ch_versions.mix(DIPLOID_SAMPLE_REHEADER.out.versions)

    emit:
        small_vcf = INDELS_SAMPLE_REHEADER.out.vcf
        small_vcfgz = INDELS_SAMPLE_REHEADER.out.vcfgz
        small_tbi = INDELS_SAMPLE_REHEADER.out.tbi
        small_csi = INDELS_SAMPLE_REHEADER.out.csi

        candidate_vcf = SV_SAMPLE_REHEADER.out.vcf
        candidate_vcfgz = SV_SAMPLE_REHEADER.out.vcfgz
        candidate_tbi = SV_SAMPLE_REHEADER.out.tbi
        candidate_csi = SV_SAMPLE_REHEADER.out.csi

        diploid_vcf = DIPLOID_SAMPLE_REHEADER.out.vcf
        diploid_vcfgz = DIPLOID_SAMPLE_REHEADER.out.vcfgz
        diploid_tbi = DIPLOID_SAMPLE_REHEADER.out.tbi
        diploid_csi = DIPLOID_SAMPLE_REHEADER.out.csi
        versions = ch_versions
}
