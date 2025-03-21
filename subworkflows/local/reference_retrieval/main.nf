/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BIODBCORE_REFSEQ      } from '../../../modules/local/biodbcore/refseq/main'

include { GUNZIP                } from '../../../modules/nf-core/gunzip/main'
include { BWA_INDEX             } from '../../../modules/nf-core/bwa/index/main'
include { TABIX_BGZIP           } from '../../../modules/nf-core/tabix/bgzip/main'
include { SAMTOOLS_FAIDX        } from '../../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_BGZIP_FAIDX  } from '../../../modules/nf-core/samtools/faidx/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow REFERENCE_RETRIEVAL {

    take:
        taxonomy_id
        outdir
        reference_genome

    main:

        BIODBCORE_REFSEQ(
            taxonomy_id,
            outdir,
            reference_genome
        )

        biodbcore_json_result = BIODBCORE_REFSEQ.out.json
            .splitJson()
            .collect()
            .map { jsonList ->
                def jsonMap = jsonList.collectEntries { [it.key, it.value] }
                return tuple(jsonMap['genome_size'] as Integer, jsonMap['genome_size_ungapped'] as Integer, file(jsonMap['genome_file']))
            }

        reference_genome_ungapped_size = biodbcore_json_result.map { it[1] }


        GUNZIP(
            BIODBCORE_REFSEQ.out.reference_genome.map { file -> tuple([id: file.simpleName.replaceFirst(/\.gz$/, '')], file) }
        )

        BWA_INDEX(
            GUNZIP.out.gunzip
        )

        TABIX_BGZIP(
            GUNZIP.out.gunzip
        )

        SAMTOOLS_BGZIP_FAIDX(
            TABIX_BGZIP.out.output,
            [[], []]
        )

        SAMTOOLS_FAIDX(
            GUNZIP.out.gunzip,
            [[], []]
        )


    emit:
        reference_genome                        = BIODBCORE_REFSEQ.out.reference_genome
        reference_genome_ungapped_size          = reference_genome_ungapped_size
        reference_genome_unzipped               = GUNZIP.out.gunzip
        reference_genome_bgzipped               = TABIX_BGZIP.out.output
        reference_genome_bwa_index              = BWA_INDEX.out.index
        reference_genome_bgzipped_index         = TABIX_BGZIP.out.gzi
        reference_genome_bgzipped_faidx         = SAMTOOLS_BGZIP_FAIDX.out.fai
        reference_genome_faidx                  = SAMTOOLS_FAIDX.out.fai
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
