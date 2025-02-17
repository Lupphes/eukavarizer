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

workflow EUKAVARIZER {

    take:
        ch_genome_file // Channel with reference genome from SEQRETRIEVAL
        ch_fastq_files // Channel with FASTQ files from SEQRETRIEVAL

    main:
        //
        // Initialize channels for versions and MultiQC
        //
        ch_versions = Channel.empty()
        ch_multiqc_files = Channel.empty()
        ch_fastq_files = ch_fastq_files
            .flatMap { fastq_list -> fastq_list.collect { fastq -> tuple([ id: fastq.baseName ], file(fastq)) } }
        ch_fastq_files.collect().view()

        //
        // MODULE: Run FastQC on Sequencing Data
        //
        FASTQC (
            ch_fastq_files
        )

        // Collect FastQC outputs for MultiQC
        ch_fastqc_report = FASTQC.out.html
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect { it[1] })
        ch_versions = ch_versions.mix(FASTQC.out.versions)

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

        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList(),
            [],
            []
        )

    emit:
        fastqc_report = ch_fastqc_report.toList() // Path to FastQC reports
        multiqc_report = MULTIQC.out.report.toList() // Path to MultiQC report
        versions = ch_versions // Path to versions.yml
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
