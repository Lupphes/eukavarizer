/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC_MULTIQC_ANALYSIS as PRE_FASTQC_MULTIQC_ANALYSIS    } from '../../../subworkflows/local/fastqc_multiqc_analysis/main'
include { FASTQC_MULTIQC_ANALYSIS as AFTER_FASTQC_MULTIQC_ANALYSIS  } from '../../../subworkflows/local/fastqc_multiqc_analysis/main'

include { FASTP         } from '../../../modules/nf-core/fastp/main'
include { BBMAP_BBDUK   } from '../../../modules/nf-core/bbmap/bbduk/main'
include { SEQTK_SAMPLE  } from '../../../modules/nf-core/seqtk/sample/main'

include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../../nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../..//nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../../local/utils_nfcore_eukavarizer_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow QUALITY_CONTROL {

    take:
        fastq_files

    main:

        PRE_FASTQC_MULTIQC_ANALYSIS(
            fastq_files
        )

        // FASTP
        FASTP(
            fastq_files,
            [],
            false,
            false,
            false
        )

        BBMAP_BBDUK(
            FASTP.out.reads
                .filter { it ->
                    if (it[1] instanceof List) {
                        // Paired-end: check both files
                        it[1][0].toFile().length() > 0 && it[1][1].toFile().length() > 0
                    } else {
                        // Single-end: check single file
                        it[1].toFile().length() > 0
                    }
                },
            []
        )

        SEQTK_SAMPLE(
            BBMAP_BBDUK.out.reads.map { meta, fastq -> tuple(meta, fastq, 0.1) }
        )

        AFTER_FASTQC_MULTIQC_ANALYSIS(
            SEQTK_SAMPLE.out.reads
        )

    emit:
        fastqc_report   = params.fastqc_flag ? AFTER_FASTQC_MULTIQC_ANALYSIS.out.fastqc_report.toList() : null      // Path to FastQC reports, or null if FastQC is not run
        multiqc_report  = params.multiqc_flag ? AFTER_FASTQC_MULTIQC_ANALYSIS.out.multiqc_report.toList() : null    // Path to MultiQC report, or null if MultiQC is not run
        processed_files = fastq_files                                                                               // Processed FASTQ files
        versions        = SEQTK_SAMPLE.out.reads                                                                    // Path to versions.yml (this is always generated)

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
