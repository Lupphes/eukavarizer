include { TABIX_BGZIP as TABIX_BGZIP_FASTQ1                     } from '../../../modules/nf-core/tabix/bgzip/main'
include { TABIX_BGZIP as TABIX_BGZIP_FASTQ2                     } from '../../../modules/nf-core/tabix/bgzip/main'

include { SEQKIT_SIZE                                           } from '../../../modules/local/seqkit/size/main'

include { SRATOOLS_FASTERQDUMP                                  } from '../../../modules/nf-core/sratools/fasterqdump/main'
include { DORADO as DORADO_FAST5                                } from '../../../modules/local/dorado/main'
include { DORADO as DORADO_POD5                                 } from '../../../modules/local/dorado/main'
include { BAM_CONVERT_SAMTOOLS as BAM_SAMTOOLS_COLLATEFASTQ     } from '../../../subworkflows/sarek/bam_convert_samtools/main'
include { BAM_CONVERT_SAMTOOLS as CRAM_SAMTOOLS_COLLATEFASTQ    } from '../../../subworkflows/sarek/bam_convert_samtools/main'

workflow SEQUENCE_FASTQ_CONVERTOR {
    take:
        samplesheet
        reference_genome_unzipped
        reference_genome_bwa_index

    main:
    input_sample_type = samplesheet.branch{
            fastq_gz:           it[0].data_type == "fastq_gz"
            fastq:              it[0].data_type == "fastq"
            bam:                it[0].data_type == "bam"
            cram:               it[0].data_type == "cram"
            sra:                it[0].data_type == "sra"
            fast5:              it[0].data_type == "fast5"
            pod5:               it[0].data_type == "pod5"
            bax_h5:             it[0].data_type == "bax_h5"
        }

        // Zip the uncompressed ones
        TABIX_BGZIP_FASTQ1(
            input_sample_type.fastq.map {meta, fastq1, _fastq2 ->
                tuple(meta, fastq1)
            }
        )

        TABIX_BGZIP_FASTQ2(
            input_sample_type.fastq.map {meta, _fastq1, fastq2 ->
                tuple(meta, fastq2)
            }
        )

        zipped_joined_ch = TABIX_BGZIP_FASTQ1.out.output.join(TABIX_BGZIP_FASTQ2.out.output, by: 0, failOnDuplicate: true, failOnMismatch: true)

        // Convert BAM or CRAM to FASTQ
        interleave_input = false  // Currently don't allow interleaved input
        BAM_SAMTOOLS_COLLATEFASTQ(
            input_sample_type.bam,
            reference_genome_unzipped.collect(),
            reference_genome_bwa_index.collect(),
            interleave_input
        )

        CRAM_SAMTOOLS_COLLATEFASTQ(
            input_sample_type.cram,
            reference_genome_unzipped.collect(),
            reference_genome_bwa_index.collect(),
            interleave_input
        )

        SRATOOLS_FASTERQDUMP(
            input_sample_type.sra,
            [],
            []
        )

        sra_file = SRATOOLS_FASTERQDUMP.out.reads.map { meta, sra_fastq ->
            tuple(meta + [single_end: sra_fastq.size() == 1], sra_fastq)
        }

        DORADO_FAST5(
            input_sample_type.fast5,
        )

        DORADO_POD5(
            input_sample_type.pod5,
        )

        fastq_gz = input_sample_type.fastq_gz.map { meta, files -> addReadgroupToMeta(meta, files) }
        zipped_joined_ch = input_sample_type.fastq.map { meta, files -> addReadgroupToMeta(meta, files) }
        sra = sra_file.map { meta, files -> addReadgroupToMeta(meta, files) }
        fast5 = DORADO_FAST5.out.fastq.map { meta, files -> addReadgroupToMeta(meta, files) }
        pod5 = DORADO_POD5.out.fastq.map { meta, files -> addReadgroupToMeta(meta, files) }

        collected_fastqs = fastq_gz
            .mix(zipped_joined_ch)
            // .mix(sra)
            .mix(fast5)
            .mix(pod5)
            .mix(BAM_SAMTOOLS_COLLATEFASTQ.out.reads)
            .mix(CRAM_SAMTOOLS_COLLATEFASTQ.out.reads)

        SEQKIT_SIZE(
            collected_fastqs
        )

        // Add median_bp to the metadata
        tagged_collected_fastqs = SEQKIT_SIZE.out.median_bp.map { meta, fastq, median_bp ->
            def length = median_bp.text.trim().toInteger()
            tuple(meta + [median_bp: length], fastq)
        }

    emit:
        tagged_collected_fastqs
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ADAPTED FUNCTIONS FROM NF-CORE/SAREK
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Add readgroup to meta and remove lane
def addReadgroupToMeta(meta, files) {
    def CN = params.seq_center ? "CN:${params.seq_center}\\t" : ''
    def flowcell = flowcellLaneFromFastq(files[0])

    // Check if flowcell ID matches
    if ( flowcell && flowcell != flowcellLaneFromFastq(files[1]) ){
        error("Flowcell ID does not match for paired reads of sample ${meta.id} - ${files}")
    }

    // If we cannot read the flowcell ID from the fastq file, then we don't use it
    def sample_lane_id = flowcell ? "${meta.flowcell}.${meta.sample}.${meta.lane}" : "${meta.sample}.${meta.lane}"

    // Don't use a random element for ID, it breaks resuming
    def read_group = "\"@RG\\tID:${sample_lane_id}\\t${CN}PU:${meta.lane}\\tSM:${meta.patient}_${meta.sample}\\tLB:${meta.sample}\\tPL:${meta.platform}\""
    meta  = meta - meta.subMap('lane') + [read_group: read_group.toString()]
    return [ meta, files ]
}

// Parse first line of a FASTQ file, return the flowcell id and lane number.
def flowcellLaneFromFastq(path) {
    // First line of FASTQ file contains sequence identifier plus optional description
    def firstLine = readFirstLineOfFastq(path)
    def flowcell_id = null

    if (!firstLine) {
        log.warn "FASTQ file(${path}): Empty or unreadable, cannot extract flowcell ID"
        return flowcell_id
    }

    // Expected format from ILLUMINA
    // cf https://en.wikipedia.org/wiki/FASTQ_format#Illumina_sequence_identifiers
    // Five fields:
    // @<instrument>:<lane>:<tile>:<x-pos>:<y-pos>...
    // Seven fields or more (from CASAVA 1.8+):
    // "@<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>..."

    def fields = firstLine ? firstLine.split(':') : []
    if (fields.size() == 5) {
        // Get the instrument name as flowcell ID
        flowcell_id = fields[0].substring(1)
    } else if (fields.size() >= 7) {
        // Get the actual flowcell ID
        flowcell_id = fields[2]
    } else if (fields.size() == 1 && fields[0].startsWith('@SRR')) {
        log.info "FASTQ file(${path}): SRA-style header detected, no flowcell ID present"
    } else if (firstLine ==~ /^@m\d+_\d+_\d+/) {
        // PacBio movie name: @m64011_191205_205309/136540764/ccs
        flowcell_id = firstLine.split('/')[0].substring(1)
        log.info "FASTQ file(${path}): Detected PacBio movie ID as flowcell: ${flowcell_id}"
    } else if (firstLine.startsWith('@read')) {
        // ONT/Dorado: No flowcell in FASTQ
        log.info "FASTQ file(${path}): Detected ONT read ID format, no flowcell ID present"
    } else if (fields.size() != 0) {
        log.warn "FASTQ file(${path}): Cannot extract flowcell ID from ${firstLine}"
    }
    return flowcell_id
}

// Get first line of a FASTQ file
def readFirstLineOfFastq(path) {
    def line = null
    try {
        InputStream stream = path.name.endsWith(".gz") ?
            new java.util.zip.GZIPInputStream(path.newInputStream()) :
            path.newInputStream()

        def reader = new BufferedReader(new InputStreamReader(stream, 'ASCII'))
        line = reader.readLine()
        if (line == null || !line.startsWith('@')) {
            log.warn "FASTQ file(${path}): Invalid or missing header line"
            return null
        }
    } catch (Exception e) {
        log.warn "FASTQ file(${path}): Error streaming"
        log.warn "${e.message}"
    }
    return line
}
