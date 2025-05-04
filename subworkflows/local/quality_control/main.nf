/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    QUALITY_CONTROL WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    This workflow performs quality control on raw FASTQ files:
    1. **PRE_FASTQC_MULTIQC_ANALYSIS** – Runs initial FastQC and MultiQC analysis.
    2. **FASTP / FASTPLONG** – Trims and filters short reads using FASTP,
        and trims/filters long reads using FASTPLONG.
    3. **BBMAP_BBDUK** – Removes adapter and contaminant sequences.
    4. **SEQTK_SAMPLE** – Subsamples reads to reduce data size for testing (if enabled).
    5. **AFTER_FASTQC_MULTIQC_ANALYSIS** – Runs final FastQC and MultiQC analysis.

    Outputs:
    - `fastq_filtered` – Processed FASTQ files.
    - `fastqc_report` – FastQC report (if enabled).
    - `multiqc_report` – MultiQC report (if enabled).
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC_MULTIQC_ANALYSIS as PRE_FASTQC_MULTIQC_ANALYSIS    } from '../../../subworkflows/local/fastqc_multiqc_analysis/main'
include { FASTQC_MULTIQC_ANALYSIS as AFTER_FASTQC_MULTIQC_ANALYSIS  } from '../../../subworkflows/local/fastqc_multiqc_analysis/main'

include { FASTP                     } from '../../../modules/nf-core/fastp/main'
include { FASTPLONG                 } from '../../../modules/local/fastplong/main'
include { BBMAP_BBDUK               } from '../../../modules/nf-core/bbmap/bbduk/main'
include { SEQTK_SAMPLE              } from '../../../modules/nf-core/seqtk/sample/main'

include { paramsSummaryMap          } from 'plugin/nf-schema'
include { paramsSummaryMultiqc      } from '../../nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML    } from '../..//nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText    } from '../../local/utils_nfcore_eukavarizer_pipeline'

workflow QUALITY_CONTROL {

    take:
        fastq_files

    main:

        PRE_FASTQC_MULTIQC_ANALYSIS(
            fastq_files
        )

        if (params.fastp_flag) {
            // FASTP
            FASTP(
                fastq_files
                // Don't use on long reads
                .filter { meta, _fastq ->
                    (
                        (meta.platform && meta.platform == 'illumina') ||
                        (!meta.platform && meta.median_bp <= params.long_read_threshold)
                    )
                },
                [],
                false,
                false,
                false
            )

            FASTPLONG(
                fastq_files
                // Don't use on short reads
                .filter { meta, _fastq ->
                    (
                        (meta.platform && meta.platform == 'ont' || meta.platform == 'pacbio') ||
                        (!meta.platform && meta.median_bp > params.long_read_threshold)
                    )
                },
                [],
                false,
                false
            )

            // FASTP Short reads & FASTPLONG Long reads
            fastp_combined_result = FASTP.out.reads.mix(FASTPLONG.out.reads)
        }
        else {
            fastp_combined_result = fastq_files
        }

        non_null_fastp_combined = fastp_combined_result.filter { it ->
            if (it[1] instanceof List) {
                // Paired-end: check both files
                it[1][0].toFile().length() > 0 && it[1][1].toFile().length() > 0
            } else {
                // Single-end: check single file
                it[1].toFile().length() > 0
            }
        }

        if (params.bbmap_bbduk_flag) {
            bbduk_result = BBMAP_BBDUK(
                non_null_fastp_combined
                    // Don't use on long reads
                    .filter { meta, _fastq ->
                        (
                            (meta.platform && meta.platform == 'illumina') ||
                            (!meta.platform && meta.median_bp <= params.long_read_threshold)
                        )
                    },
                []
            ).reads
        }
        else {
            bbduk_result = non_null_fastp_combined
        }


        // Trimmed Illumina reads and filtered long reads
        bbduk_combined_result = bbduk_result.mix(
            fastp_combined_result.filter { meta, _fastq ->
            (
                (meta.platform && meta.platform == 'ont' || meta.platform == 'pacbio') ||
                (!meta.platform && meta.median_bp > params.long_read_threshold)
            )
        })

        if (params.seqtk_flag) {
            seqtk_combined_result = SEQTK_SAMPLE(
                bbduk_combined_result.map { meta, fastq -> tuple(meta, fastq, params.seqtk_size) }
            ).reads
        } else{
            seqtk_combined_result = bbduk_combined_result
        }

        AFTER_FASTQC_MULTIQC_ANALYSIS(
            seqtk_combined_result
        )

    emit:
        fastqc_report   = params.fastqc_flag ? AFTER_FASTQC_MULTIQC_ANALYSIS.out.fastqc_report.toList() : null      // Path to FastQC reports, or null if FastQC is not run
        multiqc_report  = params.multiqc_flag ? AFTER_FASTQC_MULTIQC_ANALYSIS.out.multiqc_report.toList() : null    // Path to MultiQC report, or null if MultiQC is not run
        fastq_filtered  = seqtk_combined_result                                                                    // Processed FASTQ files
        versions        = params.multiqc_flag ? AFTER_FASTQC_MULTIQC_ANALYSIS.out.versions : null                   // Path to versions.yml (this is always generated)

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
