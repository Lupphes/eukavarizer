/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TIDDIT_SV WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    This workflow detects structural variants (SVs) using TIDDIT and processes the output:
    1. **TIDDIT_SV** – Detects SVs from BAM files using a coverage-based approach.
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
        name_tiddit = "tiddit"

        TIDDIT_BWA(
            bam_inputs
            // Only keep paired-end reads as TIDDIT does not support unpaired reads or long reads
            .filter { meta, _bam, _bai ->
                !meta.single_end
            }
            // Use BWA for short reads
            .filter { meta, _bam, _bai ->
                !params.minimap2_flag || (
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
            // Use minimap2 for long reads
            .filter { meta, _bam, _bai ->
                params.minimap2_flag && (
                    (meta.platform && meta.platform == 'illumina') ||
                    (!meta.platform && meta.median_bp <= params.long_read_threshold)
                )
            }.map { meta, bam, bai ->
                tuple(meta + [id: "${meta.id}"], bam, bai)
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
                    Channel.value(file("${projectDir}/assets/svync/${name_tiddit}.yaml"))
                )
        )

        SAMPLE_REHEADER(
            SVYNC.out.vcf,
            name_tiddit,
            ""
        )

    emit:
        vcf = SAMPLE_REHEADER.out.vcf
        vcfgz =  SAMPLE_REHEADER.out.vcfgz
        tbi = SAMPLE_REHEADER.out.tbi
        csi = SAMPLE_REHEADER.out.csi
}
