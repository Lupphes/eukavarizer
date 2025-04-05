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

include { MANTA_GERMLINE                                  } from '../../../modules/nf-core/manta/germline/main'

include { SAMPLE_REHEADER as SMALL_SAMPLE_REHEADER        } from '../../../modules/local/sample_regen/main.nf'
include { SAMPLE_REHEADER as CANDIDATE_SAMPLE_REHEADER    } from '../../../modules/local/sample_regen/main.nf'
include { SAMPLE_REHEADER as DIPLOID_SAMPLE_REHEADER      } from '../../../modules/local/sample_regen/main.nf'

include { SVYNC as SMALL_SVYNC                            } from '../../../modules/nf-core/svync/main'
include { SVYNC as CANDIDATE_SVYNC                        } from '../../../modules/nf-core/svync/main'
include { SVYNC as DIPLOID_SVYNC                          } from '../../../modules/nf-core/svync/main'

include { GUNZIP as SMALL_GUNZIP                          } from '../../../modules/nf-core/gunzip/main'
include { GUNZIP as CANDIDATE_GUNZIP                      } from '../../../modules/nf-core/gunzip/main'
include { GUNZIP as DIPLOID_GUNZIP                        } from '../../../modules/nf-core/gunzip/main'

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
                !params.minimap2_flag || meta.median_bp <= params.long_read_threshold
            }
            .map { meta, bam, bai ->
                tuple(meta + [id: "${meta.id}_${name_manta}"], bam, bai, [], [])
            },
            reference_genome_unzipped.map { meta, fasta ->
            tuple(meta + [id: "${meta.id}"], fasta)
            },
            reference_genome_faidx.map { meta, fai ->
                tuple(meta + [id: "${meta.id}"], fai)
            },
            []
        )

        SMALL_SAMPLE_REHEADER(
            MANTA_GERMLINE.out.candidate_small_indels_vcf,
            MANTA_GERMLINE.out.candidate_small_indels_vcf.map { meta, _vcf -> "${meta.id}_${name_manta_small}" },
            false
        )

        CANDIDATE_SAMPLE_REHEADER(
            MANTA_GERMLINE.out.candidate_sv_vcf,
            MANTA_GERMLINE.out.candidate_sv_vcf.map { meta, _vcf -> "${meta.id}_${name_manta_candidate}" },
            false
        )

        DIPLOID_SAMPLE_REHEADER(
            MANTA_GERMLINE.out.diploid_sv_vcf,
            MANTA_GERMLINE.out.diploid_sv_vcf.map { meta, _vcf -> "${meta.id}_${name_manta_diploid}" },
            false
        )


        SMALL_SVYNC(
            SMALL_SAMPLE_REHEADER.out.vcf
                .join(SMALL_SAMPLE_REHEADER.out.tbi, by: 0)
                .map { meta, vcf, tbi ->
                    tuple([id: "${meta.id}_svync"], vcf, tbi)
                }
                .combine(
                    Channel.value(file("${projectDir}/assets/svync/${name_manta}.yaml"))
                )
        )

        CANDIDATE_SVYNC(
            CANDIDATE_SAMPLE_REHEADER.out.vcf
                .join(CANDIDATE_SAMPLE_REHEADER.out.tbi, by: 0)
                .map { meta, vcf, tbi ->
                    tuple([id: "${meta.id}_svync"], vcf, tbi)
                }
                .combine(
                    Channel.value(file("${projectDir}/assets/svync/${name_manta}.yaml"))
                )
        )

        DIPLOID_SVYNC(
            DIPLOID_SAMPLE_REHEADER.out.vcf
                .join(DIPLOID_SAMPLE_REHEADER.out.tbi, by: 0)
                .map { meta, vcf, tbi ->
                    tuple([id: "${meta.id}_svync"], vcf, tbi)
                }
                .combine(
                    Channel.value(file("${projectDir}/assets/svync/${name_manta}.yaml"))
                )
        )

        SMALL_GUNZIP(
            SMALL_SVYNC.out.vcf
        )

        CANDIDATE_GUNZIP(
            CANDIDATE_SVYNC.out.vcf
        )

        DIPLOID_GUNZIP(
            DIPLOID_SVYNC.out.vcf
        )

    emit:
        small_vcf = SMALL_SAMPLE_REHEADER.out.vcf
        small_vcfgz = SMALL_SAMPLE_REHEADER.out.vcfgz
        small_tbi = SMALL_SAMPLE_REHEADER.out.tbi
        small_csi = SMALL_SAMPLE_REHEADER.out.csi
        svync_small_vcf = SMALL_GUNZIP.out.gunzip
        svync_small_vcfgz = SMALL_SVYNC.out.vcf
        svync_small_tbi = SMALL_SVYNC.out.tbi

        candidate_vcf = CANDIDATE_SAMPLE_REHEADER.out.vcf
        candidate_vcfgz = CANDIDATE_SAMPLE_REHEADER.out.vcfgz
        candidate_tbi = CANDIDATE_SAMPLE_REHEADER.out.tbi
        candidate_csi = CANDIDATE_SAMPLE_REHEADER.out.csi
        svync_candidate_vcf = CANDIDATE_GUNZIP.out.gunzip
        svync_candidate_vcfgz = CANDIDATE_SVYNC.out.vcf
        svync_candidate_tbi = CANDIDATE_SVYNC.out.tbi

        diploid_vcf = DIPLOID_SAMPLE_REHEADER.out.vcf
        diploid_vcfgz = DIPLOID_SAMPLE_REHEADER.out.vcfgz
        diploid_tbi = DIPLOID_SAMPLE_REHEADER.out.tbi
        diploid_csi = DIPLOID_SAMPLE_REHEADER.out.csi
        svync_diploid_vcf = DIPLOID_GUNZIP.out.gunzip
        svync_diploid_vcfgz = DIPLOID_SVYNC.out.vcf
        svync_diploid_tbi = DIPLOID_SVYNC.out.tbi
}
