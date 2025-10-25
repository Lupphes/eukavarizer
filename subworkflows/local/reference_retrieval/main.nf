/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: REFERENCE_RETRIEVAL
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Downloads, processes, and indexes reference genome files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Description:
        This subworkflow retrieves reference genomes from RefSeq or uses provided files,
        then prepares them for analysis by creating necessary indexes for alignment and
        variant calling tools.

    Processing Steps:
        1. Download reference genome from RefSeq (BIODBCORE_REFSEQ)
        2. Decompress reference genome (GUNZIP)
        3. Generate BWA or BWA-MEM2 alignment index
        4. Generate Minimap2 alignment index (if enabled)
        5. Compress genome with bgzip (TABIX_BGZIP)
        6. Create FASTA indexes for both compressed and uncompressed versions (SAMTOOLS_FAIDX)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Input Channels (take):
        taxonomy_id               val(id)                NCBI taxonomy identifier
        outdir                    path(dir)              Output directory path
        reference_genome          path(fasta)            Reference genome file (optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Output Channels (emit):
        reference_genome          path(fasta.gz)         Original reference genome
        reference_genome_unzipped tuple(meta, fasta)     Uncompressed FASTA file
        reference_genome_bgzipped tuple(meta, fasta.gz)  bgzipped FASTA file
        reference_genome_bwa_index tuple(meta, index)    BWA/BWA-MEM2 index files
        reference_genome_minimap_index tuple(meta, mmi)  Minimap2 index file
        reference_genome_faidx    tuple(meta, fai)       FASTA index for uncompressed
        reference_genome_bgzipped_faidx tuple(meta, fai) FASTA index for bgzipped
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Author:   Ondřej Sloup (Lupphes)
    Contact:  ondrej.sloup@protonmail.com
    GitHub:   @Lupphes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BIODBCORE_REFSEQ      } from '../../../modules/local/biodbcore/refseq/main'

include { GUNZIP                } from '../../../modules/nf-core/gunzip/main'
include { BWAMEM2_INDEX         } from '../../../modules/nf-core/bwamem2/index/main'
include { BWA_INDEX             } from '../../../modules/nf-core/bwa/index/main'
include { MINIMAP2_INDEX        } from '../../../modules/nf-core/minimap2/index/main'
include { TABIX_BGZIP           } from '../../../modules/nf-core/tabix/bgzip/main'
include { SAMTOOLS_FAIDX        } from '../../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_BGZIP_FAIDX  } from '../../../modules/nf-core/samtools/faidx/main'

workflow REFERENCE_RETRIEVAL {

    take:
        taxonomy_id
        outdir
        reference_genome

    main:
        ch_versions = channel.empty()

        BIODBCORE_REFSEQ(
            taxonomy_id,
            outdir,
            reference_genome
        )

        biodbcore_json_result = BIODBCORE_REFSEQ.out.json
            .splitJson()
            .collect()
            .map { jsonList ->
                def jsonMap = jsonList.collectEntries { entry -> [entry.key, entry.value] }
                return tuple(jsonMap['genome_size'] as Integer, jsonMap['genome_size_ungapped'] as Integer, file(jsonMap['reference_genome']))
            }

        reference_genome_ungapped_size = biodbcore_json_result.map { result -> result[1] }

        if (reference_genome != []) {
            // User provided reference
            // In production: reference_genome is a channel from channel.fromPath()
            // In tests: reference_genome is a file() object
            // Solution: Use flatMap to handle both - it unwraps channels and works with single values
            reference_genome_input = channel.from(reference_genome)
                .flatMap { x -> (x instanceof Collection) ? x : [x] }
                .map { file -> tuple([id: file.simpleName.replaceFirst(/\.gz$/, '')], file) }
        } else {
            // BIODBCORE output - channel that needs collect/flatten/map
            reference_genome_input = BIODBCORE_REFSEQ.out.reference_genome
                .collect()
                .flatten()
                .map { file -> tuple([id: file.simpleName.replaceFirst(/\.gz$/, '')], file) }
        }

        GUNZIP(
            reference_genome_input
        )

        bwa_index = (params.bwamem2 ? BWAMEM2_INDEX(GUNZIP.out.gunzip) : BWA_INDEX(GUNZIP.out.gunzip))

        if (params.minimap2_flag) {
            minimap_index_ch = MINIMAP2_INDEX(
                GUNZIP.out.gunzip
            ).index
        }

        TABIX_BGZIP(
            GUNZIP.out.gunzip
        )

        SAMTOOLS_BGZIP_FAIDX(
            TABIX_BGZIP.out.output,
            [[], []],
            false
        )

        SAMTOOLS_FAIDX(
            GUNZIP.out.gunzip,
            [[], []],
            false
        )


        ch_versions = ch_versions.mix(BIODBCORE_REFSEQ.out.versions)
        ch_versions = ch_versions.mix(GUNZIP.out.versions)
        ch_versions = ch_versions.mix(bwa_index.versions)
        if (params.minimap2_flag) {
            ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)
        }
        ch_versions = ch_versions.mix(TABIX_BGZIP.out.versions)
        ch_versions = ch_versions.mix(SAMTOOLS_BGZIP_FAIDX.out.versions)
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    emit:
        reference_genome                        = reference_genome_input
        reference_genome_ungapped_size          = reference_genome_ungapped_size
        reference_genome_unzipped               = GUNZIP.out.gunzip
        reference_genome_bgzipped               = TABIX_BGZIP.out.output
        reference_genome_bwa_index              = bwa_index.index
        reference_genome_minimap_index          = params.minimap2_flag ? minimap_index_ch : channel.empty()
        reference_genome_bgzipped_index         = TABIX_BGZIP.out.gzi
        reference_genome_bgzipped_faidx         = SAMTOOLS_BGZIP_FAIDX.out.fai
        reference_genome_faidx                  = SAMTOOLS_FAIDX.out.fai
        versions                                = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
