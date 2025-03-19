/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC as AFTER_FASTQC    } from '../modules/nf-core/fastqc/main'
include { MULTIQC as AFTER_MULTIQC  } from '../modules/nf-core/multiqc/main'
include { FASTQC as PRE_FASTQC      } from '../modules/nf-core/fastqc/main'
include { MULTIQC as PRE_MULTIQC  } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_eukavarizer_pipeline'

include { FASTP } from '../modules/nf-core/fastp/main'
include { BBMAP_BBDUK } from '../modules/nf-core/bbmap/bbduk/main'
include { SEQTK_SAMPLE } from '../modules/nf-core/seqtk/sample/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ANALYSIS_FAST_MULTI {

    take:
        ch_fastq_files  // Channel with FASTQ files from SEQRETRIEVAL
        debug_flag      // Debug flag
        fastqc_flag     // FastQC flag
        multiqc_flag    // MultiQC flag

    main:
        //
        // Initialize channels for versions and MultiQC
        //
        ch_versions = Channel.empty()
        ch_multiqc_files = Channel.empty()

        if (fastqc_flag) {
            combined_fastq_files = ch_fastq_files
                .map { fileList ->
                    if (!(fileList instanceof List) || fileList.size() != 3) {
                        log.error("Invalid FASTQ input format: ${fileList}")
                        return null
                    }

                    def sample_id = fileList[0].toString()
                    def paired_end_files = fileList[1] ?: []
                    def single_end_files = fileList[2] ?: []

                    def files
                    def is_paired

                    if (paired_end_files.size() == 2) {
                        files = paired_end_files
                        is_paired = true
                    } else if (single_end_files.size() == 1) {
                        files = single_end_files
                        is_paired = false
                    } else {
                        log.warn("Invalid file count for sample ${sample_id}: Paired -> ${paired_end_files.size()}, Single -> ${single_end_files.size()}")
                        return null
                    }

                    def meta = [id: sample_id, single_end: !is_paired]
                    return tuple(meta, files)
                }
                .filter { it != null }





            //
            // MODULE: Run FastQC on Sequencing Data
            //
            PRE_FASTQC(
                combined_fastq_files
            )

            // Collect FastQC outputs for MultiQC
            ch_fastqc_report = PRE_FASTQC.out.html
            ch_multiqc_files = ch_multiqc_files.mix(PRE_FASTQC.out.zip.collect { it[1] })
            ch_versions = ch_versions.mix(PRE_FASTQC.out.versions)

        }

        //
        // Collate and save software versions
        //
        softwareVersionsToYAML(ch_versions)
            .collectFile(
                storeDir: "${params.outdir}/pipeline_info",
                name: 'nf_core_eukavarizer_software_mqc_versions.yml',
                sort: true,
                newLine: true
            ).set { ch_collated_versions }


        if (multiqc_flag) {
            //
            // MODULE: MultiQC
            //
            ch_multiqc_config = Channel.fromPath(
                "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
            ch_multiqc_custom_config = params.multiqc_config ?
                Channel.fromPath(params.multiqc_config, checkIfExists: true) :
                Channel.empty()
            ch_multiqc_logo = params.multiqc_logo ?
                Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
                Channel.empty()

            summary_params = paramsSummaryMap(
                workflow, parameters_schema: "nextflow_schema.json")
            ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
            ch_multiqc_files = ch_multiqc_files.mix(
                ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))

            ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
                file(params.multiqc_methods_description, checkIfExists: true) :
                file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
            ch_methods_description = Channel.value(
                methodsDescriptionText(ch_multiqc_custom_methods_description))

            ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
            ch_multiqc_files = ch_multiqc_files.mix(
                ch_methods_description.collectFile(
                    name: 'methods_description_mqc.yaml',
                    sort: true
                )
            )

            PRE_MULTIQC(
                ch_multiqc_files.collect(),
                ch_multiqc_config.toList(),
                ch_multiqc_custom_config.toList(),
                ch_multiqc_logo.toList(),
                [], // Placeholder if you need to add more config files
                []  // Placeholder for any additional parameters
            )

            if (debug_flag) {
                ch_multiqc_files.collect().view { "DEBUG: MultiQC files -> ${it}" }
            }
        }


                // Normalize FASTQ input into a structured format
        ch_formatted_fastq = ch_fastq_files
            .map { fileList ->
                if (!(fileList instanceof List) || fileList.size() != 3) {
                    log.error("Invalid FASTQ input format: ${fileList}")
                    return null
                }

                def sample_id = fileList[0].toString()
                def paired_end_files = fileList[1] ?: []
                def single_end_files = fileList[2] ?: []

                def files = []
                def is_paired = false

                if (paired_end_files.size() == 2) {
                    files = paired_end_files
                    is_paired = true
                } else if (single_end_files.size() == 1) {
                    files = single_end_files
                    is_paired = false
                } else {
                    log.warn("Invalid file count for sample ${sample_id}: Paired -> ${paired_end_files.size()}, Single -> ${single_end_files.size()}")
                    return null
                }

                def meta = [ id: sample_id, single_end: !is_paired ]
                return tuple(meta, files)
            }
            .filter { it != null }


        // FASTP
        fastp = FASTP(
            ch_formatted_fastq,
            [],
            false,
            false,
            false
        )

        fastp.reads.collect().view { "DEBUG: FASTP reads -> ${it}" }
        // fastp.reads.collect().view { "DEBUG: ETOO -> ${it}" }
        // fastp.reads_fail.collect().view { "DEBUG: FASTP reads_fail -> ${it}" }
        // fastp.reads_merged.collect().view { "DEBUG: FASTP reads_merged -> ${it}" }

        filtered = fastp.reads
            .filter { it ->
                if (it[1] instanceof List) {
                    // Paired-end: check both files
                    it[1][0].toFile().length() > 0 && it[1][1].toFile().length() > 0
                } else {
                    // Single-end: check single file
                    it[1].toFile().length() > 0
                }
            }


        filtered.collect().view { "DEBUG: FASTP filtered reads -> ${it}" }

        // BBMAP_BBDUK

        bbmap = BBMAP_BBDUK(
            filtered,
            []
        )

        bbmap.reads.collect().view { "DEBUG: BBMAP_BBDUK reads -> ${it}" }

        input = bbmap.reads.map { meta, fastq -> tuple(meta, fastq, 0.1) }

        input.collect().view { "DEBUG: SEQTK_SAMPLE input -> ${it}" }

        // SEQTK

        seqtk = SEQTK_SAMPLE(
            input
        )

        seqtk.reads.collect().view { "DEBUG: SEQTK_SAMPLE reads -> ${it}" }


        ch_versions_after = Channel.empty()
        ch_multiqc_files_after = Channel.empty()


        // FASTQC (After SEQTK)
        if (fastqc_flag) {

        combined_fastq_files_after = seqtk.reads
            .map { inputo ->
                if (!(inputo instanceof List) || inputo.size() != 2) {
                    log.error("Invalid FASTQ inputo format: ${inputo}")
                    return null
                }

                def meta1 = inputo[0]         // Meta info (id, single_end)
                def reads2 = inputo[1]        // Reads (single or list of paired)

                // Ensure reads2 are properly formatted
                if (meta1.single_end && reads2 instanceof Path) {
                    // Single-end case
                    return tuple(meta1, [reads2])
                } else if (!meta1.single_end && reads2 instanceof List && reads2.size() == 2) {
                    // Paired-end case
                    return tuple(meta1, reads2)
                } else {
                    log.warn("Invalid file count for sample ${meta1.id}: Reads -> ${reads2}")
                    return null
                }
            }
            .filter { it != null }


            AFTER_FASTQC(
                combined_fastq_files_after
            )

            // Collect FastQC outputs for MultiQC
            ch_fastqc_report_after = AFTER_FASTQC.out.html
            ch_multiqc_files_after = ch_multiqc_files_after.mix(AFTER_FASTQC.out.zip.collect { it[1] })
            ch_versions_after = ch_versions_after.mix(AFTER_FASTQC.out.versions)


            if (debug_flag) {
                seqtk.reads.collect().view { "DEBUG: SEQTK_SAMPLE reads -> ${it}" }
            }
        }
        softwareVersionsToYAML(ch_versions_after)
            .collectFile(
                storeDir: "${params.outdir}/pipeline_info",
                name: 'nf_core_eukavarizer_software_mqc_versions.yml',
                sort: true,
                newLine: true
            ).set { ch_collated_versions_after }

        // MultiQC (After AFTER_FASTQC)
        if (multiqc_flag) {
            ch_multiqc_config_after = Channel.fromPath(
                "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
            ch_multiqc_custom_config_after = params.multiqc_config ?
                Channel.fromPath(params.multiqc_config, checkIfExists: true) :
                Channel.empty()
            ch_multiqc_logo_after = params.multiqc_logo ?
                Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
                Channel.empty()

            summary_params_after = paramsSummaryMap(
                workflow, parameters_schema: "nextflow_schema.json")
            ch_workflow_summary_after = Channel.value(paramsSummaryMultiqc(summary_params_after))
            ch_multiqc_files_after = ch_multiqc_files_after.mix(
                ch_workflow_summary_after.collectFile(name: 'workflow_summary_mqc.yaml'))

            ch_multiqc_custom_methods_description_after = params.multiqc_methods_description ?
                file(params.multiqc_methods_description, checkIfExists: true) :
                file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
            ch_methods_description_after = Channel.value(
                methodsDescriptionText(ch_multiqc_custom_methods_description_after))

            ch_multiqc_files_after = ch_multiqc_files_after.mix(ch_collated_versions_after)
            ch_multiqc_files_after = ch_multiqc_files_after.mix(
                ch_methods_description_after.collectFile(
                    name: 'methods_description_mqc.yaml',
                    sort: true
                )
            )

            AFTER_MULTIQC(
                ch_multiqc_files_after.collect(),
                ch_multiqc_config_after.toList(),
                ch_multiqc_custom_config_after.toList(),
                ch_multiqc_logo_after.toList(),
                [], // Placeholder if you need to add more config files
                []  // Placeholder for any additional parameters
            )

            if (debug_flag) {
                ch_multiqc_files_after.collect().view { "DEBUG: MultiQC files -> ${it}" }
            }
        }

        seqtk.reads.collect().view { "DEBUG: SEQTK_SAMPLE reads -> ${it}" }
        ch_fastq_files.collect().view { "DEBUG: FASTQ files -> ${it}" }

        formatted_reads = seqtk.reads.collect { it ->
            if (it[1] instanceof List) {
                // Paired-end reads
                [it[0].id, it[1], []]
            } else {
                // Single-end reads
                [it[0].id, [], [it[1]]]
            }
        }


    emit:
        fastqc_report = fastqc_flag ? ch_fastqc_report_after.toList() : null // Path to FastQC reports, or null if FastQC is not run
        multiqc_report = multiqc_flag ? AFTER_MULTIQC.out.report.toList() : null // Path to MultiQC report, or null if MultiQC is not run
        processed_files = ch_fastq_files // Processed FASTQ files
        versions = formatted_reads                              // Path to versions.yml (this is always generated)

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
