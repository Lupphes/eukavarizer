/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: QUALITY_CONTROL
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Performs comprehensive quality control, trimming, and filtering of FASTQ files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Description:
        This subworkflow applies quality control measures to raw sequencing reads,
        including adapter trimming, quality filtering, contamination removal, and
        optional subsampling. It generates QC reports before and after processing.

    Processing Steps:
        1. Run initial FastQC and MultiQC analysis (PRE_FASTQC_MULTIQC_ANALYSIS)
        2. Trim and filter short reads with FASTP or long reads with FASTPLONG
        3. Remove adapter and contaminant sequences (BBMAP_BBDUK)
        4. Subsample reads for testing if enabled (SEQTK_SAMPLE)
        5. Run final FastQC and MultiQC analysis (AFTER_FASTQC_MULTIQC_ANALYSIS)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Input Channels (take):
        fastq_files               tuple(meta, fastq)     Input FASTQ files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Output Channels (emit):
        fastq_filtered            tuple(meta, fastq)     Filtered and trimmed FASTQ files
        multiqc_report            path(html)             MultiQC quality control report
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Author:   Ondřej Sloup (Lupphes)
    Contact:  ondrej.sloup@protonmail.com
    GitHub:   @Lupphes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC_MULTIQC_ANALYSIS as PRE_FASTQC_MULTIQC_ANALYSIS    } from '../../../subworkflows/local/fastqc_multiqc_analysis/main'
include { FASTQC_MULTIQC_ANALYSIS as AFTER_FASTQC_MULTIQC_ANALYSIS  } from '../../../subworkflows/local/fastqc_multiqc_analysis/main'

include { FASTP                     } from '../../../modules/nf-core/fastp/main'
include { FASTPLONG                 } from '../../../modules/nf-core/fastplong/main'
include { BBMAP_BBDUK               } from '../../../modules/nf-core/bbmap/bbduk/main'
include { SEQTK_SAMPLE              } from '../../../modules/nf-core/seqtk/sample/main'

workflow QUALITY_CONTROL {

    take:
        fastq_files

    main:
        ch_versions = Channel.empty()

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
                }
                .map { meta, reads -> tuple(meta, reads, []) },
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
                }
                .map { meta, reads -> tuple(meta, reads, []) },
                false,
                false,
                false
            )

            // FASTP Short reads & FASTPLONG Long reads
            fastp_combined_result = FASTP.out.reads.mix(FASTPLONG.out.reads)
        }
        else {
            fastp_combined_result = fastq_files
        }

        // Filter out samples with no reads (empty files)
        // This is critical to prevent downstream errors and ensure accurate sample counting
        non_null_fastp_combined = fastp_combined_result.filter { it ->
            def meta = it[0]
            def files = it[1]
            def hasReads = false

            if (files instanceof List) {
                // Paired-end: check both files
                hasReads = files[0].toFile().length() > 0 && files[1].toFile().length() > 0
            } else {
                // Single-end: check single file
                hasReads = files.toFile().length() > 0
            }

            if (!hasReads) {
                log.warn "Sample ${meta.id} (${meta.sample}) has 0 reads after quality control and will be excluded from downstream analysis (alignment, merging, SV calling)"
            }

            return hasReads
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

        ch_versions = ch_versions.mix(PRE_FASTQC_MULTIQC_ANALYSIS.out.versions)
        ch_versions = ch_versions.mix(AFTER_FASTQC_MULTIQC_ANALYSIS.out.versions)
        if (params.fastp_flag) {
            ch_versions = ch_versions.mix(FASTP.out.versions.first())
            ch_versions = ch_versions.mix(FASTPLONG.out.versions.first())
        }
        if (params.bbmap_bbduk_flag) {
            ch_versions = ch_versions.mix(BBMAP_BBDUK.out.versions.first())
        }
        if (params.seqtk_flag) {
            ch_versions = ch_versions.mix(SEQTK_SAMPLE.out.versions.first())
        }

    emit:
        fastqc_report   = params.fastqc_flag ? AFTER_FASTQC_MULTIQC_ANALYSIS.out.fastqc_report.toList() : null      // Path to FastQC reports, or null if FastQC is not run
        multiqc_report  = params.multiqc_flag ? AFTER_FASTQC_MULTIQC_ANALYSIS.out.multiqc_report.toList() : null    // Path to MultiQC report, or null if MultiQC is not run
        fastq_filtered  = seqtk_combined_result                                                                    // Processed FASTQ files
        versions        = ch_versions                                                                              // Software versions

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
