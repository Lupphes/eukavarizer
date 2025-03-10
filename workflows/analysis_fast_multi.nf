/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_eukavarizer_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ANALYSIS_FAST_MULTI {

    take:
        ch_fastq_files  // Channel with FASTQ files from SEQRETRIEVAL
        debug_flag   // Debug flag
        fastqc_flag  // FastQC flag
        multiqc_flag // MultiQC flag

    main:
        //
        // Initialize channels for versions and MultiQC
        //
        ch_versions = Channel.empty()
        ch_multiqc_files = Channel.empty()

        // Normalize FASTQ input into a structured format
        ch_formatted_fastq = ch_fastq_files
            .map { fileList ->
                def files = fileList instanceof List ? fileList.flatten() : [fileList] // Ensure list format
                def sample_id = files[0].toString().tokenize('/')[-1].replaceAll(/(_\d+)?\.fastq\.gz$/, "")
                def is_paired = files.size() == 2
                return tuple([ id: sample_id, single_end: !is_paired ], files)
            }

        if (debug_flag) {
            ch_formatted_fastq.view { "DEBUG: Formatted FASTQ -> ${it}" }
        }

        if (fastqc_flag) {
            //
            // MODULE: Run FastQC on Sequencing Data
            //
            FASTQC(
                ch_formatted_fastq
            )

            // Collect FastQC outputs for MultiQC
            ch_fastqc_report = FASTQC.out.html
            ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect { it[1] })
            ch_versions = ch_versions.mix(FASTQC.out.versions)

            if (debug_flag) {
                ch_formatted_fastq.view { "DEBUG: Formatted FASTQ -> ${it}" }
            }
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

            MULTIQC(
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


    emit:
        fastqc_report = fastqc_flag ? ch_fastqc_report.toList() : null // Path to FastQC reports, or null if FastQC is not run
        multiqc_report = multiqc_flag ? MULTIQC.out.report.toList() : null // Path to MultiQC report, or null if MultiQC is not run
        versions = ch_versions // Path to versions.yml (this is always generated)

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
