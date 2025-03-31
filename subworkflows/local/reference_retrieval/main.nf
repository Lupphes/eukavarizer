/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    REFERENCE_RETRIEVAL WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    This workflow retrieves and processes reference genomes:
    1. **BIODBCORE_REFSEQ** – Downloads reference genome from RefSeq.
    2. **GUNZIP** – Decompresses the reference genome.
    3. **BWAMEM2_INDEX** or **BWA_INDEX** – Creates BWA2 index for the reference genome.
    4. **MINIMAP2_INDEX** – Creates MINIMAP2 index for the reference genome.
    5. **TABIX_BGZIP** – Compresses the reference genome with bgzip.
    6. **SAMTOOLS_FAIDX** – Creates a FASTA index for the uncompressed genome.
    7. **SAMTOOLS_BGZIP_FAIDX** – Creates a FASTA index for the bgzipped genome.

    Outputs:
    - `reference_genome` – Original reference genome file.
    - `reference_genome_unzipped` – Unzipped reference genome.
    - `reference_genome_bgzipped` – Bgzipped reference genome.
    - `reference_genome_bwa_index` – BWA index.
    - `reference_genome_faidx` – FASTA index for the uncompressed genome.
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
                return tuple(jsonMap['genome_size'] as Integer, jsonMap['genome_size_ungapped'] as Integer, file(jsonMap['reference_genome']))
            }

        reference_genome_ungapped_size = biodbcore_json_result.map { it[1] }
        reference_genome_input = (reference_genome != [] ? reference_genome : BIODBCORE_REFSEQ.out.reference_genome).collect().flatten()

        GUNZIP(
            reference_genome_input.map { file -> tuple([id: file.simpleName.replaceFirst(/\.gz$/, '')], file) }
        )

        bwa_index = (params.bwamem2 ? BWAMEM2_INDEX(GUNZIP.out.gunzip) : BWA_INDEX(GUNZIP.out.gunzip))

        MINIMAP2_INDEX(
            GUNZIP.out.gunzip
        )

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


    emit:
        reference_genome                        = reference_genome_input
        reference_genome_ungapped_size          = reference_genome_ungapped_size
        reference_genome_unzipped               = GUNZIP.out.gunzip
        reference_genome_bgzipped               = TABIX_BGZIP.out.output
        reference_genome_bwa_index              = bwa_index.index
        reference_genome_minimap_index          = MINIMAP2_INDEX.out.index
        reference_genome_bgzipped_index         = TABIX_BGZIP.out.gzi
        reference_genome_bgzipped_faidx         = SAMTOOLS_BGZIP_FAIDX.out.fai
        reference_genome_faidx                  = SAMTOOLS_FAIDX.out.fai
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
