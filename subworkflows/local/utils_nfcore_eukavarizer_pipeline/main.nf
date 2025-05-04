/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PIPELINE INITIALISATION AND COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PIPELINE_INITIALISATION:
    - Sets up the pipeline by validating input parameters and configuration.
    - Prints the pipeline version and configuration information.
    - Optionally parses the samplesheet and validates input formats.

    PIPELINE_COMPLETION:
    - Handles the pipeline's completion tasks.
    - Sends completion emails (HTML/plaintext) and generates summary logs.
    - Sends IM notifications if configured.
    - Catches and logs any pipeline errors.

    Outputs:
    - `versions` – Software version summary (currently empty channel placeholder).
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BIODBCORE_ENA         } from '../../../modules/local/biodbcore/ena/main'

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'


workflow PIPELINE_INITIALISATION {

    take:
    version             // boolean: Display version and exit
    validate_params     // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs     // boolean: Do not use coloured log outputs
    nextflow_cli_args   // array: List of positional nextflow CLI args
    outdir              // string: The output directory where the results will be saved
    input               // string: Path to input samplesheet


    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Create channel from input file provided through params.input
    //

    if (input) {
        // Adapted from nf-core/sarek
        Channel
            .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
            .map { meta, fastq_1, fastq_2, bam, cram, sra, bax_h5, fast5, pod5 ->
                [ meta.patient + meta.sample, [meta, fastq_1, fastq_2, bam, cram, sra, bax_h5, fast5, pod5] ]
            }.tap { ch_with_patient_sample } // save the channel
            .groupTuple() //group by patient_sample to get all lanes
            .map { patient_sample, ch_items ->
                [ patient_sample, ch_items.size() ]
            }.combine(ch_with_patient_sample, by: 0) // for each entry add numLanes
            .map { patient_sample, num_lanes, ch_items ->
                def (meta, fastq_1, fastq_2, bam, cram, sra, bax_h5, fast5, pod5) = ch_items
                def isZipped = { file -> file?.name.endsWith('.gz') }

                if (meta.lane && fastq_1 && fastq_2) {
                    def zipped1 = isZipped(fastq_1)
                    def zipped2 = isZipped(fastq_2)

                    if (zipped1 != zipped2) {
                        throw new IllegalArgumentException("FASTQ pair has inconsistent compression: ${fastq_1.name} and ${fastq_2.name}")
                    }
                    def zipped = zipped1 && zipped2

                    meta = meta + [id: "${meta.sample}-${meta.lane}".toString(), single_end: false, data_type : zipped ? "fastq_gz" : "fastq", num_lanes: num_lanes.toInteger(), size: 1, platform: meta.platform]
                    return [ meta, [ fastq_1, fastq_2 ] ]

                } else if (fastq_1 && fastq_2) {
                    def zipped1 = isZipped(fastq_1)
                    def zipped2 = isZipped(fastq_2)

                    if (zipped1 != zipped2) {
                        throw new IllegalArgumentException("FASTQ pair has inconsistent compression: ${fastq_1.name} and ${fastq_2.name}")
                    }
                    def zipped = zipped1 && zipped2

                    meta = meta + [id: meta.sample, single_end: false, data_type : zipped ? "fastq_gz" : "fastq", platform: meta.platform]
                    return [ meta - meta.subMap('lane'), [ fastq_1, fastq_2 ] ]

                } else if (meta.lane && fastq_1 && !fastq_2) {
                    def zipped = isZipped(fastq_1)
                    meta = meta + [id: "${meta.sample}-${meta.lane}".toString(), single_end: true, data_type : zipped ? "fastq_gz" : "fastq", num_lanes: num_lanes.toInteger(), size: 1, platform: meta.platform]
                    return [ meta, [ fastq_1 ] ]

                } else if (fastq_1 && !fastq_2) {
                    def zipped = isZipped(fastq_1)
                    meta = meta + [id: meta.sample, single_end: true, data_type : zipped ? "fastq_gz" : "fastq", platform: meta.platform]
                    return [ meta - meta.subMap('lane'), [ fastq_1 ] ]

                } else if (meta.lane && bam) {
                    meta            = meta + [id: "${meta.sample}-${meta.lane}".toString()]
                    def CN          = params.seq_center ? "CN:${params.seq_center}\\t" : ''
                    def read_group  = "\"@RG\\tID:${meta.sample}_${meta.lane}\\t${CN}PU:${meta.lane}\\tSM:${meta.patient}_${meta.sample}\\tLB:${meta.sample}\\tPL:${meta.platform}\""

                    meta            = meta - meta.subMap('lane') + [read_group: read_group.toString(), single_end: true, data_type: 'bam', num_lanes: num_lanes.toInteger(), size: 1, platform: meta.platform]
                    return [ meta - meta.subMap('lane'), bam ]

                } else if (bam) {
                    meta = meta + [id: meta.sample, data_type: 'bam', platform: meta.platform]

                    return [ meta - meta.subMap('lane'), bam ]

                } else if (meta.lane && cram) {
                    meta            = meta + [id: "${meta.sample}-${meta.lane}".toString()]
                    def CN          = params.seq_center ? "CN:${params.seq_center}\\t" : ''
                    def read_group  = "\"@RG\\tID:${meta.sample}_${meta.lane}\\t${CN}PU:${meta.lane}\\tSM:${meta.patient}_${meta.sample}\\tLB:${meta.sample}\\tPL:${meta.platform}\""

                    meta            = meta - meta.subMap('lane') + [read_group: read_group.toString(), single_end: true, data_type: 'cram', num_lanes: num_lanes.toInteger(), size: 1, platform: meta.platform]
                    return [ meta - meta.subMap('lane'), cram ]

                } else if (cram) {
                    meta = meta + [id: meta.sample, data_type: 'cram', platform: meta.platform]
                    return [ meta - meta.subMap('lane'), cram ]

                } else if (meta.lane && sra) {
                    // Single-end SRA will be parsed to fastq and then updated
                    meta = meta + [id: "${meta.sample}-${meta.lane}".toString(), single_end: true, data_type : "sra", num_lanes: num_lanes.toInteger(), size: 1, platform: meta.platform]
                    return [ meta, sra ]

                } else if (sra) {
                    // Single-end SRA will be parsed to fastq and then updated
                    meta = meta + [id: meta.sample, single_end: true, data_type : "sra", platform: meta.platform]
                    return [ meta - meta.subMap('lane'), sra ]

                } else if (meta.lane && fast5) {
                    meta = meta - meta.subMap('lane') + [id: "${meta.sample}-${meta.lane}".toString(), single_end: true, data_type: 'fast5', num_lanes: num_lanes.toInteger(), size: 1, platform: meta.platform]
                    return [ meta - meta.subMap('lane'), fast5 ]

                } else if (fast5) {
                    meta = meta + [id: meta.sample, single_end: true, data_type: 'fast5', platform: meta.platform]
                    return [ meta - meta.subMap('lane'), fast5 ]

                } else if (meta.lane && bax_h5) {
                    meta = meta - meta.subMap('lane') + [id: "${meta.sample}-${meta.lane}".toString(), single_end: true, data_type: 'bax_h5', num_lanes: num_lanes.toInteger(), size: 1, platform: meta.platform]
                    return [ meta - meta.subMap('lane'), bax_h5 ]

                } else if (bax_h5) {
                    meta = meta + [id: meta.sample, single_end: true, data_type: 'bax_h5', platform: meta.platform]
                    return [ meta - meta.subMap('lane'), bax_h5 ]

                } else {
                    error("Samplesheet does not contain a sequence. Please check your samplesheet or adjust it accordingly.")
                }
        }.set { ch_samplesheet }
    } else {
        log.info "No samplesheet provided — will search for input FASTQs automatically"
        Channel.empty().set { [] }

        //TODO: This should just return samplesheet
        BIODBCORE_ENA(
            params.taxonomy_id,
            params.reference_genome,
            outdir,
            params.library_strategy.map     { it.join(' ').trim() },
            params.instrument_platform.map  { it.join(' ').trim() },
            params.minimum_coverage,
            params.maximum_coverage,
            params.max_results,
            params.assembly_quality,
        )
    }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}

//
// Generate methods description for MultiQC
//
def toolCitationText() {
    def citation_text = [
        "Tools used in the workflow included:",
        "FastQC (Andrews 2010),",
        "MultiQC (Ewels et al. 2016),",
        "BWA (Li and Durbin 2009),",
        "BWA-MEM2 (Vasimuddin et al. 2019),",
        "Minimap2 (Li 2018),",
        "Samtools (Li et al. 2009),",
        "BCFtools (Danecek et al. 2021),",
        "Fastp (Chen et al. 2018),",
        "fastplong (Chen 2023),",
        "BBDuk (Bushnell 2014),",
        "Seqtk (Li 2012),",
        "DELLY (Rausch et al. 2012),",
        "GRIDSS (Cameron et al. 2017),",
        "TIDDIT (Eisfeldt et al. 2017),",
        "DYSGU (Cawte et al. 2022),",
        "Sniffles (Sedlazeck et al. 2018),",
        "CuteSV (Jiang et al. 2020),",
        "SVABA (Wala et al. 2018),",
        "Manta (Chen et al. 2016),",
        "SURVIVOR (Jeffares et al. 2017),",
        "SVYnc (Nieuwkoop et al. 2023),",
        "Varify (Sloup et al. 2024),",
        "SRATools (Leinonen et al. 2011),",
        "BioDbCore (Sloup 2023),",
        "Coreutils/Gunzip (GNU Project)"
    ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    def reference_text = [
        "<li>Andrews S. (2010) FastQC. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/</li>",
        "<li>Ewels P. et al. (2016) MultiQC. Bioinformatics, 32(19), 3047–3048. https://doi.org/10.1093/bioinformatics/btw354</li>",
        "<li>Li H., Durbin R. (2009) BWA. Bioinformatics, 25(14), 1754–1760. https://doi.org/10.1093/bioinformatics/btp324</li>",
        "<li>Vasimuddin M. et al. (2019) BWA-MEM2. arXiv:1907.12931</li>",
        "<li>Li H. (2018) Minimap2. Bioinformatics, 34(18), 3094–3100. https://doi.org/10.1093/bioinformatics/bty191</li>",
        "<li>Li H. et al. (2009) SAMtools. Bioinformatics, 25(16), 2078–2079. https://doi.org/10.1093/bioinformatics/btp352</li>",
        "<li>Danecek P. et al. (2021) BCFtools. GigaScience, 10(2), giab008. https://doi.org/10.1093/gigascience/giab008</li>",
        "<li>Chen S., Zhou Y., Chen Y., Gu J. (2018) fastp. Bioinformatics, 34(17), i884–i890. https://doi.org/10.1093/bioinformatics/bty560</li>",
        "<li>Chen S. (2023) Ultrafast one-pass FASTQ data preprocessing using fastp. iMeta 2: e107. https://doi.org/10.1002/imt2.107</li>",
        "<li>Bushnell B. (2014) BBMap. https://sourceforge.net/projects/bbmap/</li>",
        "<li>Li H. (2012) seqtk. https://github.com/lh3/seqtk</li>",
        "<li>Rausch T. et al. (2012) DELLY. Bioinformatics, 28(18), i333–i339. https://doi.org/10.1093/bioinformatics/bts378</li>",
        "<li>Cameron DL. et al. (2017) GRIDSS. Genome Research, 27(12), 2050–2060. https://doi.org/10.1101/gr.222109.117</li>",
        "<li>Eisfeldt J. et al. (2017) TIDDIT. F1000Research, 6:664. https://doi.org/10.12688/f1000research.11501.1</li>",
        "<li>Cawte N. et al. (2022) Dysgu. Bioinformatics, 38(9), 2465–2471. https://doi.org/10.1093/bioinformatics/btac144</li>",
        "<li>Sedlazeck FJ. et al. (2018) Sniffles. Nature Methods, 15, 461–468. https://doi.org/10.1038/s41592-018-0001-7</li>",
        "<li>Jiang T. et al. (2020) CuteSV. Genome Biology, 21:189. https://doi.org/10.1186/s13059-020-02107-y</li>",
        "<li>Wala JA. et al. (2018) SvABA. Genome Research, 28(4), 581–591. https://doi.org/10.1101/gr.221028.117</li>",
        "<li>Chen X. et al. (2016) Manta. Bioinformatics, 32(8), 1220–1222. https://doi.org/10.1093/bioinformatics/btv710</li>",
        "<li>Jeffares DC. et al. (2017) SURVIVOR. Nat Commun, 8, 14061. https://doi.org/10.1038/ncomms14061</li>",
        "<li>Nicolas Vannieuwkerke (2024) SVYNC. https://github.com/nvnieuwk/svync</li>",
        "<li>Sloup O. et al. (2024) Varify. https://github.com/Lupphes/Varify</li>",
        "<li>Leinonen R. et al. (2011) The Sequence Read Archive. Nucleic Acids Res., 39(Database issue): D19–D21. https://doi.org/10.1093/nar/gkq1019</li>",
        "<li>Sloup O. (2023) BioDbCore: ENA/RefSeq retrieval module. https://github.com/luppo/biodbcore</li>",
        "<li>GNU Project. coreutils (including gunzip). https://www.gnu.org/software/coreutils/</li>"
    ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()

    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

