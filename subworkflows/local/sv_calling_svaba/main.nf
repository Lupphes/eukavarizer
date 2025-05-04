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
                tuple(meta + [id: "${meta.id}"], bam, bai, [], [])
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
                    Channel.value(file("${projectDir}/assets/svync/${name_svaba}.yaml"))
                )
        )

        SAMPLE_REHEADER(
            SVYNC.out.vcf,
            name_svaba,
            ""
        )

    emit:
        vcf = SAMPLE_REHEADER.out.vcf
        vcfgz =  SAMPLE_REHEADER.out.vcfgz
        tbi = SAMPLE_REHEADER.out.tbi
        csi = SAMPLE_REHEADER.out.csi
}
