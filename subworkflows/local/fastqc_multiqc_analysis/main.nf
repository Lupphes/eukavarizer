/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FASTQC_MULTIQC_ANALYSIS WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    This workflow runs quality control and summarises the results:
    1. **FASTQC** – Runs FastQC on raw FASTQ files (if enabled).
    2. **MULTIQC** – Aggregates reports from FastQC and other sources (if enabled).

    Outputs:
    - `fastqc_report` – FastQC report (if enabled).
    - `multiqc_report` – MultiQC report (if enabled).
    - `versions` – Collected software versions for reproducibility.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC        } from '../../../modules/nf-core/fastqc/main'
include { MULTIQC       } from '../../../modules/nf-core/multiqc/main'

include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../../nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../..//nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../../local/utils_nfcore_eukavarizer_pipeline'

workflow FASTQC_MULTIQC_ANALYSIS {

    take:
        fastq_files

    main:
        ch_versions         = Channel.empty()
        ch_multiqc_files    = Channel.empty()

        if (params.fastqc_flag) {
            //
            // MODULE: Run FastQC on Sequencing Data
            //
            FASTQC(
                fastq_files
            )

            // Collect FastQC outputs for MultiQC
            ch_fastqc_report = FASTQC.out.html
            ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect { it[1] })
            ch_versions = ch_versions.mix(FASTQC.out.versions)
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


        if (params.multiqc_flag) {
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
                [],
                []
            )
        }

    emit:
        fastqc_report = params.fastqc_flag ? ch_fastqc_report.toList() : null       // Path to FastQC reports, or null if FastQC is not run
        multiqc_report = params.multiqc_flag ? MULTIQC.out.report.toList() : null   // Path to MultiQC report, or null if MultiQC is not run
        versions = ch_versions                                                      // Processed FASTQ files
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
