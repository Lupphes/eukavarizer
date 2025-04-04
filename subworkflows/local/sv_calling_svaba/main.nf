/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SVABA WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    This workflow detects structural variants (SVs) and indels using SVABA:
    1. **SVABA** – Calls SVs and indels using local assembly.
    2. **SAMPLE_REHEADER** – Renames and reheaders the VCF output.
    3. **SVYNC** – Synchronizes and refines SV calls.
    4. **GUNZIP** – Decompresses the final VCF file.

    Outputs:
    - `vcf`           – Reheaded VCF file.
    - `vcfgz`         – Gzipped reheaded VCF file.
    - `tbi`           – Tabix index (.tbi) for the reheaded VCF.
    - `csi`           – CSI index (.csi) for the reheaded VCF (if generated).
    - `svync_vcf`     – Decompressed synchronized VCF file.
    - `svync_vcfgz`   – Gzipped synchronized VCF file.
    - `svync_tbi`     – Tabix index for the synchronized VCF.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SVABA             } from '../../../modules/nf-core/svaba/main'
include { SAMPLE_REHEADER   } from '../../../modules/local/sample_regen/main.nf'
include { SVYNC             } from '../../../modules/nf-core/svync/main'
include { GUNZIP            } from '../../../modules/nf-core/gunzip/main'

workflow SV_CALLING_SVABA {
    take:
        bam_inputs
        reference_genome_unzipped
        reference_genome_faidx
        reference_genome_bwa_index

    main:
        name_svaba = "svaba"

        SVABA(
            bam_inputs.filter { meta, _bam, _bai ->
                !params.minimap2_flag || meta.median_bp <= params.long_read_threshold
            }.map { meta, bam, bai ->
                tuple(meta + [id: "${meta.id}_svaba"], bam, bai, [], [])
            },
            reference_genome_unzipped,
            reference_genome_faidx,
            reference_genome_bwa_index,
            [[],[]],
            [[],[]],
            [[],[]],
        )

        SAMPLE_REHEADER(
            SVABA.out.sv,
            SVABA.out.sv.map { meta, _vcf -> "${meta.id}_${name_svaba}" },
            true
        )

        SVYNC(
            SAMPLE_REHEADER.out.vcf
                .join(SAMPLE_REHEADER.out.tbi, by: 0)
                .map { meta, vcf, tbi ->
                    tuple([id: "${meta.id}_svync"], vcf, tbi)
                }
                .combine(
                    Channel.value(file("${projectDir}/assets/svync/${name_svaba}.yaml"))
                )
        )

        GUNZIP(
            SVYNC.out.vcf
        )

    emit:
        vcf = SAMPLE_REHEADER.out.vcf
        vcfgz = SAMPLE_REHEADER.out.vcfgz
        tbi = SAMPLE_REHEADER.out.tbi
        csi = SAMPLE_REHEADER.out.csi
        svync_vcf = GUNZIP.out.gunzip
        svync_vcfgz = SVYNC.out.vcf
        svync_tbi = SVYNC.out.tbi
}
