/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MANTA_GERMLINE WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    This workflow calls structural variants (SVs) using MANTA for germline samples:
    1. **MANTA_GERMLINE** – Detects large deletions, inversions, and duplications.
    2. **SAMPLE_REHEADER** – Reheaders and renames output VCF files.
    3. **SVYNC** – Synchronizes and refines SV calls.
    4. **GUNZIP** – Decompresses the final VCF files.

    Outputs:
    - `small_vcf`              – Reheaded VCF file for small indels.
    - `small_vcfgz`            – Gzipped VCF for small indels.
    - `small_tbi`              – Tabix index for small VCF.
    - `small_csi`              – CSI index for small VCF (if generated).
    - `svync_small_vcf`        – Decompressed synchronized small VCF.
    - `svync_small_vcfgz`      – Gzipped synchronized small VCF.
    - `svync_small_tbi`        – Tabix index for synchronized small VCF.

    - `candidate_vcf`          – Reheaded VCF file for candidate SVs.
    - `candidate_vcfgz`        – Gzipped VCF for candidate SVs.
    - `candidate_tbi`          – Tabix index for candidate VCF.
    - `candidate_csi`          – CSI index for candidate VCF (if generated).
    - `svync_candidate_vcf`    – Decompressed synchronized candidate VCF.
    - `svync_candidate_vcfgz`  – Gzipped synchronized candidate VCF.
    - `svync_candidate_tbi`    – Tabix index for synchronized candidate VCF.

    - `diploid_vcf`            – Reheaded VCF file for diploid SVs.
    - `diploid_vcfgz`          – Gzipped VCF for diploid SVs.
    - `diploid_tbi`            – Tabix index for diploid VCF.
    - `diploid_csi`            – CSI index for diploid VCF (if generated).
    - `svync_diploid_vcf`      – Decompressed synchronized diploid VCF.
    - `svync_diploid_vcfgz`    – Gzipped synchronized diploid VCF.
    - `svync_diploid_tbi`      – Tabix index for synchronized diploid VCF.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MANTA_GERMLINE                                } from '../../../modules/nf-core/manta/germline/main'

include { SAMPLE_REHEADER as INDELS_SAMPLE_REHEADER     } from '../../../modules/local/sample_regen/main.nf'
include { SAMPLE_REHEADER as SV_SAMPLE_REHEADER         } from '../../../modules/local/sample_regen/main.nf'
include { SAMPLE_REHEADER as DIPLOID_SAMPLE_REHEADER    } from '../../../modules/local/sample_regen/main.nf'

include { SVYNC as INDEL_SVYNC                          } from '../../../modules/nf-core/svync/main'
include { SVYNC as SV_SVYNC                             } from '../../../modules/nf-core/svync/main'
include { SVYNC as DIPLOID_SVYNC                        } from '../../../modules/nf-core/svync/main'

include { GUNZIP as INDEL_GUNZIP                        } from '../../../modules/nf-core/gunzip/main'
include { GUNZIP as SV_GUNZIP                           } from '../../../modules/nf-core/gunzip/main'
include { GUNZIP as DIPLOID_GUNZIP                      } from '../../../modules/nf-core/gunzip/main'

workflow SV_CALLING_MANTA {
    take:
        bam_inputs
        reference_genome_unzipped
        reference_genome_faidx

    main:
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

        // TODO: Make proper fix fox manta when VCF file completelly empty and don't rely first on REHEADER
        INDELS_SAMPLE_REHEADER(
            MANTA_GERMLINE.out.candidate_small_indels_vcf,
            name_manta,
            "-${name_manta_small}"
        )

        SV_SAMPLE_REHEADER(
            MANTA_GERMLINE.out.candidate_sv_vcf,
            name_manta,
            "-${name_manta_candidate}"
        )

        DIPLOID_SAMPLE_REHEADER(
            MANTA_GERMLINE.out.diploid_sv_vcf,
            name_manta,
            "-${name_manta_diploid}"
        )

        INDEL_SVYNC(
            INDELS_SAMPLE_REHEADER.out.vcf
                .join(INDELS_SAMPLE_REHEADER.out.tbi, by: 0)
                .map { meta, vcf, tbi ->
                    tuple(meta + [id: "${meta.id}-svync"], vcf, tbi)
                }
                .combine(
                    Channel.value(file("${projectDir}/assets/svync/${name_manta}.yaml"))
                )
        )

        SV_SVYNC(
            SV_SAMPLE_REHEADER.out.vcf
                .join(SV_SAMPLE_REHEADER.out.tbi, by: 0)
                .map { meta, vcf, tbi ->
                    tuple(meta + [id: "${meta.id}-svync"], vcf, tbi)
                }
                .combine(
                    Channel.value(file("${projectDir}/assets/svync/${name_manta}.yaml"))
                )
        )

        DIPLOID_SVYNC(
            DIPLOID_SAMPLE_REHEADER.out.vcf
                .join(DIPLOID_SAMPLE_REHEADER.out.tbi, by: 0)
                .map { meta, vcf, tbi ->
                    tuple(meta + [id: "${meta.id}-svync"], vcf, tbi)
                }
                .combine(
                    Channel.value(file("${projectDir}/assets/svync/${name_manta}.yaml"))
                )
        )

        INDEL_GUNZIP(
            INDEL_SVYNC.out.vcf
        )

        SV_GUNZIP(
            SV_SVYNC.out.vcf
        )

        DIPLOID_GUNZIP(
            DIPLOID_SVYNC.out.vcf
        )

    emit:
        small_vcf = INDEL_GUNZIP.out.gunzip
        small_vcfgz = INDEL_SVYNC.out.vcf
        small_tbi = INDEL_SVYNC.out.tbi
        small_csi = INDELS_SAMPLE_REHEADER.out.csi

        candidate_vcf = SV_GUNZIP.out.gunzip
        candidate_vcfgz = SV_SVYNC.out.vcf
        candidate_tbi = SV_SVYNC.out.tbi
        candidate_csi = SV_SAMPLE_REHEADER.out.csi

        diploid_vcf = DIPLOID_GUNZIP.out.gunzip
        diploid_vcfgz = DIPLOID_SVYNC.out.vcf
        diploid_tbi = DIPLOID_SVYNC.out.tbi
        diploid_csi = DIPLOID_SAMPLE_REHEADER.out.csi
}
