/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    EUKAVARIZER WORKFLOW
    Structural Variant (SV) calling pipeline for short- and long-read sequencing data.
    This workflow integrates multiple SV callers:
    - Delly, Manta, GRIDSS, Dysgu, TIDDIT, SvABA – for short reads
    - Sniffles, CuteSV – for long reads
    - Delly, Dysgu can also be used for long reads

    The outputs of this workflow include all synchronized and filtered VCFs along with
    their associated index files and metadata.

    Outputs:
    - `vcf_list`     – List of synchronized VCF files (decompressed).
    - `vcfgz_list`   – List of synchronized and compressed VCF files (.vcf.gz).
    - `tbi_list`     – List of Tabix index files (.tbi) corresponding to the compressed VCFs.
    - `meta_list`    – List of metadata maps for each variant call, used for downstream merging.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SV_CALLING_DELLY      } from '../../subworkflows/local/sv_calling_delly/main.nf'
include { SV_CALLING_MANTA      } from '../../subworkflows/local/sv_calling_manta/main.nf'
include { SV_CALLING_GRIDSS     } from '../../subworkflows/local/sv_calling_gridss/main.nf'
include { SV_CALLING_DYSGU      } from '../../subworkflows/local/sv_calling_dysgu/main.nf'
include { SV_CALLING_TIDDIT     } from '../../subworkflows/local/sv_calling_tiddit/main.nf'
include { SV_CALLING_SVABA      } from '../../subworkflows/local/sv_calling_svaba/main.nf'
include { SV_CALLING_SNIFFLES   } from '../../subworkflows/local/sv_calling_sniffles/main.nf'
include { SV_CALLING_CUTESV     } from '../../subworkflows/local/sv_calling_cutesv/main.nf'

workflow EUKAVARIZER {

    take:
        fastq_bam                           // SEQUENCE_PROCESSOR.out.bam_files
        fastq_bam_indexes                   // SEQUENCE_PROCESSOR.out.bam_indexes
        reference_genome_bgzipped           // REFERENCE_RETRIEVAL.out.reference_genome_bgzipped,
        reference_genome_faidx              // REFERENCE_RETRIEVAL.out.reference_genome_faidx
        reference_genome_bwa_index          // REFERENCE_RETRIEVAL.out.reference_genome_bwa_index
        reference_genome_bgzipped_faidx     // REFERENCE_RETRIEVAL.out.reference_genome_bgzipped_faidx
        reference_genome_unzipped           // REFERENCE_RETRIEVAL.out.reference_genome_unzipped
        reference_genome_minimap_index      // REFERENCE_RETRIEVAL.out.reference_genome_minimap_index

    main:
        bam_inputs = fastq_bam
            .join(fastq_bam_indexes, by: 0)
            .map { meta, bam, bai -> tuple(meta, bam, bai) }

        vcf_list    = Channel.value([])
        vcfgz_list  = Channel.value([])
        tbi_list    = Channel.value([])
        meta_list   = Channel.value([])

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        DELLY CALL
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        if (params.delly_flag) {
            SV_CALLING_DELLY(
                fastq_bam,
                fastq_bam_indexes,
                reference_genome_unzipped.collect(),
                reference_genome_faidx.collect()
            )

            vcf_list = vcf_list.concat(SV_CALLING_DELLY.out.svync_vcf)
            vcfgz_list = vcfgz_list.concat(SV_CALLING_DELLY.out.svync_vcfgz)
            tbi_list = tbi_list.concat(SV_CALLING_DELLY.out.svync_tbi)
            meta_list = meta_list.concat(SV_CALLING_DELLY.out.svync_vcf.map { meta, _vcf -> meta })
        }
        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        MANTA
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        if (params.manta_flag) {
            SV_CALLING_MANTA(
                bam_inputs,
                reference_genome_unzipped.collect(),
                reference_genome_faidx.collect()
            )

            vcf_list = vcf_list.concat(SV_CALLING_MANTA.out.svync_small_vcf)
            vcfgz_list = vcfgz_list.concat(SV_CALLING_MANTA.out.svync_small_vcfgz)
            tbi_list = tbi_list.concat(SV_CALLING_MANTA.out.svync_small_tbi)
            meta_list = meta_list.concat(SV_CALLING_MANTA.out.svync_small_vcf.map { meta, _vcf -> meta })

            vcf_list = vcf_list.concat(SV_CALLING_MANTA.out.svync_candidate_vcf)
            vcfgz_list = vcfgz_list.concat(SV_CALLING_MANTA.out.svync_candidate_vcfgz)
            tbi_list = tbi_list.concat(SV_CALLING_MANTA.out.svync_candidate_tbi)
            meta_list = meta_list.concat(SV_CALLING_MANTA.out.svync_candidate_vcf.map { meta, _vcf -> meta })

            vcf_list = vcf_list.concat(SV_CALLING_MANTA.out.svync_diploid_vcf)
            vcfgz_list = vcfgz_list.concat(SV_CALLING_MANTA.out.svync_diploid_vcfgz)
            tbi_list = tbi_list.concat(SV_CALLING_MANTA.out.svync_diploid_tbi)
            meta_list = meta_list.concat(SV_CALLING_MANTA.out.svync_diploid_vcf.map { meta, _vcf -> meta })
        }
        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        GRIDSS
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        if (params.gridss_flag) {
            SV_CALLING_GRIDSS(
                bam_inputs,
                reference_genome_unzipped.collect(),
                reference_genome_faidx.collect(),
                reference_genome_bwa_index.collect(),
                reference_genome_minimap_index.collect()
            )

            vcf_list = vcf_list.concat(SV_CALLING_GRIDSS.out.svync_vcf)
            vcfgz_list = vcfgz_list.concat(SV_CALLING_GRIDSS.out.svync_vcfgz)
            tbi_list = tbi_list.concat(SV_CALLING_GRIDSS.out.svync_tbi)
            meta_list = meta_list.concat(SV_CALLING_GRIDSS.out.svync_vcf.map { meta, _vcf -> meta })
        }
        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        DYSGU
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        if (params.dysgu_flag) {
            SV_CALLING_DYSGU(
                bam_inputs,
                reference_genome_unzipped.collect(),
                reference_genome_faidx.collect()
            )

            vcf_list = vcf_list.concat(SV_CALLING_DYSGU.out.svync_vcf)
            vcfgz_list = vcfgz_list.concat(SV_CALLING_DYSGU.out.svync_vcfgz)
            tbi_list = tbi_list.concat(SV_CALLING_DYSGU.out.svync_tbi)
            meta_list = meta_list.concat(SV_CALLING_DYSGU.out.svync_vcf.map { meta, _vcf -> meta })
        }
        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        TIDDIT_SV
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        if (params.tiddit_flag) {
            SV_CALLING_TIDDIT(
                bam_inputs,
                reference_genome_unzipped.collect(),
                reference_genome_bwa_index.collect(),
                reference_genome_minimap_index.collect()
            )

            vcf_list = vcf_list.concat(SV_CALLING_TIDDIT.out.svync_vcf)
            vcfgz_list = vcfgz_list.concat(SV_CALLING_TIDDIT.out.svync_vcfgz)
            tbi_list = tbi_list.concat(SV_CALLING_TIDDIT.out.svync_tbi)
            meta_list = meta_list.concat(SV_CALLING_TIDDIT.out.svync_vcf.map { meta, _vcf -> meta })
        }
        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        SVABA
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        if (params.svaba_flag && !params.bwamem2) {
            SV_CALLING_SVABA(
                bam_inputs,
                reference_genome_unzipped.collect(),
                reference_genome_faidx.collect(),
                reference_genome_bwa_index.collect()
            )

            vcf_list = vcf_list.concat(SV_CALLING_SVABA.out.svync_vcf)
            vcfgz_list = vcfgz_list.concat(SV_CALLING_SVABA.out.svync_vcfgz)
            tbi_list = tbi_list.concat(SV_CALLING_SVABA.out.svync_tbi)
            meta_list = meta_list.concat(SV_CALLING_SVABA.out.svync_vcf.map { meta, _vcf -> meta })
        }

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        SNIFFLES
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        if (params.sniffles_flag) {
            SV_CALLING_SNIFFLES(
                bam_inputs,
                reference_genome_bgzipped.collect()
            )

            vcf_list = vcf_list.concat(SV_CALLING_SNIFFLES.out.svync_vcf)
            vcfgz_list = vcfgz_list.concat(SV_CALLING_SNIFFLES.out.svync_vcfgz)
            tbi_list = tbi_list.concat(SV_CALLING_SNIFFLES.out.svync_tbi)
            meta_list = meta_list.concat(SV_CALLING_SNIFFLES.out.svync_vcf.map { meta, _vcf -> meta })
        }

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        CUTESV
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        if (params.cutesv_flag) {
            SV_CALLING_CUTESV(
                bam_inputs,
                reference_genome_unzipped.collect()
            )

            vcf_list = vcf_list.concat(SV_CALLING_CUTESV.out.svync_vcf)
            vcfgz_list = vcfgz_list.concat(SV_CALLING_CUTESV.out.svync_vcfgz)
            tbi_list = tbi_list.concat(SV_CALLING_CUTESV.out.svync_tbi)
            meta_list = meta_list.concat(SV_CALLING_CUTESV.out.svync_vcf.map { meta, _vcf -> meta })
        }


    emit:
        vcf_list = vcf_list
        vcfgz_list = vcfgz_list
        tbi_list = tbi_list
        meta_list = meta_list
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
