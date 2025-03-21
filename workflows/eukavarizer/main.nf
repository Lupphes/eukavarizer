/*
 * Pipeline: EUKAVARIZER
 * This workflow demonstrates a SV calling pipeline
 * that currently runs DELLY Manta, GRIDSS, Dysgu, TIDDIT, SVABA,
 * Sniffles, and CuteSV.
 */

include { SV_CALLING_DELLY      } from '../../subworkflows/local/sv_calling_delly/main.nf'
include { SV_CALLING_MANTA      } from '../../subworkflows/local/sv_calling_manta/main.nf'
include { SV_CALLING_GRIDSS     } from '../../subworkflows/local/sv_calling_gridss/main.nf'
include { SV_CALLING_DYSGU      } from '../../subworkflows/local/sv_calling_dysgu/main.nf'
include { SV_CALLING_TIDDIT     } from '../../subworkflows/local/sv_calling_tiddit/main.nf'
include { SV_CALLING_SVABA      } from '../../subworkflows/local/sv_calling_svaba/main.nf'
include { SV_CALLING_SNIFFLES   } from '../../subworkflows/local/sv_calling_sniffles/main.nf'
include { SV_CALLING_CUTESV     } from '../../subworkflows/local/sv_calling_cutesv/main.nf'

/*
    ──────────────────────────────────────────────────────────
    WORKFLOW DEFINITION
    ──────────────────────────────────────────────────────────
*/

workflow EUKAVARIZER {

    take:
        fastq_bam                           // SEQUENCE_PROCESSOR.out.bam_files                         [ch_bam_files]
        fastq_bam_indexes                   // SEQUENCE_PROCESSOR.out.bam_indexes                       [ch_bam_indexes]
        reference_genome_bgzipped           // REFERENCE_RETRIEVAL.out.reference_genome_bgzipped,       [ch_genome_file]
        reference_genome_faidx              // REFERENCE_RETRIEVAL.out.reference_genome_faidx           [ch_fasta_index]
        reference_genome_bwa_index          // REFERENCE_RETRIEVAL.out.reference_genome_bwa_index       [ch_bwa_index]
        reference_genome_bgzipped_faidx     // REFERENCE_RETRIEVAL.out.reference_genome_bgzipped_faidx  [ch_fasta_index_gz]
        reference_genome_unzipped           // REFERENCE_RETRIEVAL.out.reference_genome_unzipped        [ch_fasta_file]

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
                reference_genome_bgzipped,
                reference_genome_faidx
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
                reference_genome_bgzipped,
                reference_genome_bgzipped_faidx
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
                reference_genome_unzipped,
                reference_genome_faidx,
                reference_genome_bwa_index
            )

            vcf_list = vcf_list.concat(SV_CALLING_GRIDSS.out.svync_vcf)
            vcfgz_list = vcfgz_list.concat(SV_CALLING_GRIDSS.out.svync_vcfgz)
            tbi_list = tbi_list.concat(SV_CALLING_GRIDSS.out.svync_tbi)
            meta_list = meta_list.concat(SV_CALLING_GRIDSS.out.svync_vcf.map { meta, _vcf -> meta })
        }
        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        DYSGU (Not currently supported)
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        if (params.dysgu_flag) {
            SV_CALLING_DYSGU(
                bam_inputs,
                reference_genome_unzipped,
                reference_genome_faidx
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
                reference_genome_bgzipped,
                reference_genome_bwa_index
            )

            vcf_list = vcf_list.concat(SV_CALLING_TIDDIT.out.svync_vcf)
            vcfgz_list = vcfgz_list.concat(SV_CALLING_TIDDIT.out.svync_vcfgz)
            tbi_list = tbi_list.concat(SV_CALLING_TIDDIT.out.svync_tbi)
            meta_list = meta_list.concat(SV_CALLING_TIDDIT.out.svync_vcf.map { meta, _vcf -> meta })
        }
        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        SVABA (Not currently supported)
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        if (params.svaba_flag) {
            SV_CALLING_SVABA(
                bam_inputs,
                reference_genome_unzipped,
                reference_genome_faidx,
                reference_genome_bwa_index
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
                reference_genome_bgzipped
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
                reference_genome_bgzipped
            )

            vcf_list = vcf_list.concat(SV_CALLING_CUTESV.out.svync_vcf)
            vcfgz_list = vcfgz_list.concat(SV_CALLING_CUTESV.out.svync_vcfgz)
            tbi_list = tbi_list.concat(SV_CALLING_CUTESV.out.svync_tbi)
            meta_list = meta_list.concat(SV_CALLING_CUTESV.out.svync_vcf.map { meta, _vcf -> meta })
        }


    emit:
        delly_vcf                       = params.delly_flag  ?  (SV_CALLING_DELLY?.out?.vcf                     ?: null) : null
        delly_vcfgz                     = params.delly_flag  ?  (SV_CALLING_DELLY?.out?.vcfgz                   ?: null) : null
        delly_tbi                       = params.delly_flag  ?  (SV_CALLING_DELLY?.out?.tbi                     ?: null) : null
        delly_csi                       = params.delly_flag  ?  (SV_CALLING_DELLY?.out?.csi                     ?: null) : null

        delly_svync_vcf                 = params.delly_flag  ?  (SV_CALLING_DELLY?.out?.svync_vcfgz             ?: null) : null
        delly_svync_vcfgz               = params.delly_flag  ?  (SV_CALLING_DELLY?.out?.svync_vcfgz             ?: null) : null
        delly_svync_tbi                 = params.delly_flag  ?  (SV_CALLING_DELLY?.out?.svync_tbi               ?: null) : null

        manta_small_vcf                 = params.manta_flag  ?  (SV_CALLING_MANTA?.out?.small_vcf               ?: null) : null
        manta_small_vcfgz               = params.manta_flag  ?  (SV_CALLING_MANTA?.out?.small_vcfgz             ?: null) : null
        manta_small_tbi                 = params.manta_flag  ?  (SV_CALLING_MANTA?.out?.small_tbi               ?: null) : null
        manta_small_csi                 = params.manta_flag  ?  (SV_CALLING_MANTA?.out?.small_csi               ?: null) : null

        manta_svync_vcf_small           = params.manta_flag  ?  (SV_CALLING_MANTA?.out?.svync_small_vcf         ?: null) : null
        manta_svync_vcfgz_small         = params.manta_flag  ?  (SV_CALLING_MANTA?.out?.svync_small_vcfgz       ?: null) : null
        manta_svync_tbi_small           = params.manta_flag  ?  (SV_CALLING_MANTA?.out?.svync_small_tbi         ?: null) : null

        manta_candidate_vcf             = params.manta_flag  ?  (SV_CALLING_MANTA?.out?.candidate_vcf           ?: null) : null
        manta_candidate_vcfgz           = params.manta_flag  ?  (SV_CALLING_MANTA?.out?.candidate_vcfgz         ?: null) : null
        manta_candidate_tbi             = params.manta_flag  ?  (SV_CALLING_MANTA?.out?.candidate_tbi           ?: null) : null
        manta_candidate_csi             = params.manta_flag  ?  (SV_CALLING_MANTA?.out?.candidate_csi           ?: null) : null

        manta_svync_vcf_candidate       = params.manta_flag  ?  (SV_CALLING_MANTA?.out?.svync_candidate_vcf     ?: null) : null
        manta_svync_vcfgz_candidate     = params.manta_flag  ?  (SV_CALLING_MANTA?.out?.svync_candidate_vcfgz   ?: null) : null
        manta_svync_tbi_candidate       = params.manta_flag  ?  (SV_CALLING_MANTA?.out?.svync_candidate_tbi     ?: null) : null

        manta_diploid_vcf               = params.manta_flag  ?  (SV_CALLING_MANTA?.out?.diploid_vcf             ?: null) : null
        manta_diploid_vcfgz             = params.manta_flag  ?  (SV_CALLING_MANTA?.out?.diploid_vcfgz           ?: null) : null
        manta_diploid_tbi               = params.manta_flag  ?  (SV_CALLING_MANTA?.out?.diploid_tbi             ?: null) : null
        manta_diploid_csi               = params.manta_flag  ?  (SV_CALLING_MANTA?.out?.diploid_csi             ?: null) : null

        mant_svync_vcf_diploid          = params.manta_flag  ?  (SV_CALLING_MANTA?.out?.svync_diploid_vcf       ?: null) : null
        mant_svync_vcfgz_diploid        = params.manta_flag  ?  (SV_CALLING_MANTA?.out?.svync_diploid_vcfgz     ?: null) : null
        mant_svync_tbi_diploid          = params.manta_flag  ?  (SV_CALLING_MANTA?.out?.svync_diploid_tbi       ?: null) : null

        gridss_vcfs                     = params.gridss_flag  ?  (SV_CALLING_GRIDSS?.out?.vcf                   ?: null) : null
        gridss_vcfgz                    = params.gridss_flag  ?  (SV_CALLING_GRIDSS?.out?.vcfgz                 ?: null) : null
        gridss_tbi                      = params.gridss_flag  ?  (SV_CALLING_GRIDSS?.out?.tbi                   ?: null) : null
        gridss_csi                      = params.gridss_flag  ?  (SV_CALLING_GRIDSS?.out?.csi                   ?: null) : null

        gridss_svync_vcf                = params.gridss_flag ? (SV_CALLING_GRIDSS?.out?.svync_vcf               ?: null) : null
        gridss_svync_vcfgz              = params.gridss_flag ? (SV_CALLING_GRIDSS?.out?.svync_vcfgz             ?: null) : null
        gridss_svync_tbi                = params.gridss_flag ? (SV_CALLING_GRIDSS?.out?.svync_tbi               ?: null) : null

        dysgu_vcfs                      = params.dysgu_flag ? (SV_CALLING_DYSGU?.out?.vcf                       ?: null) : null
        dysgu_vcfgz                     = params.dysgu_flag ? (SV_CALLING_DYSGU?.out?.vcfgz                     ?: null) : null
        dysgu_tbi                       = params.dysgu_flag ? (SV_CALLING_DYSGU?.out?.tbi                       ?: null) : null
        dysgu_csi                       = params.dysgu_flag ? (SV_CALLING_DYSGU?.out?.csi                       ?: null) : null

        dysgu_svync_vcf                 = params.dysgu_flag ? (SV_CALLING_DYSGU?.out?.svync_vcf                 ?: null) : null
        dysgu_svync_vcfgz               = params.dysgu_flag ? (SV_CALLING_DYSGU?.out?.svync_vcfgz               ?: null) : null
        dysgu_svync_tbi                 = params.dysgu_flag ? (SV_CALLING_DYSGU?.out?.svync_tbi                 ?: null) : null

        tiddit_vcfs                     = params.tiddit_flag ? (SV_CALLING_TIDDIT?.out?.vcf                     ?: null) : null
        tiddit_vcfgz                    = params.tiddit_flag ? (SV_CALLING_TIDDIT?.out?.vcfgz                   ?: null) : null
        tiddit_tbi                      = params.tiddit_flag ? (SV_CALLING_TIDDIT?.out?.tbi                     ?: null) : null
        tiddit_csi                      = params.tiddit_flag ? (SV_CALLING_TIDDIT?.out?.csi                     ?: null) : null

        tiddit_svync_vcf                = params.tiddit_flag ? (SV_CALLING_TIDDIT?.out?.svync_vcf               ?: null) : null
        tiddit_svync_vcfgz              = params.tiddit_flag ? (SV_CALLING_TIDDIT?.out?.svync_vcfgz             ?: null) : null
        tiddit_svync_tbi                = params.tiddit_flag ? (SV_CALLING_TIDDIT?.out?.svync_tbi               ?: null) : null

        svaba_vcfs                      = params.svaba_flag ? (SV_CALLING_SVABA?.out?.vcf                       ?: null) : null
        svaba_vcfgz                     = params.svaba_flag ? (SV_CALLING_SVABA?.out?.vcfgz                     ?: null) : null
        svaba_tbi                       = params.svaba_flag ? (SV_CALLING_SVABA?.out?.tbi                       ?: null) : null
        svaba_csi                       = params.svaba_flag ? (SV_CALLING_SVABA?.out?.csi                       ?: null) : null

        svaba_svync_vcf                 = params.svaba_flag ? (SV_CALLING_SVABA?.out?.svync_vcf                 ?: null) : null
        svaba_svync_vcfgz               = params.svaba_flag ? (SV_CALLING_SVABA?.out?.svync_vcfgz               ?: null) : null
        svaba_svync_tbi                 = params.svaba_flag ? (SV_CALLING_SVABA?.out?.svync_tbi                 ?: null) : null

        sniffles_vcfs                   = params.sniffles_flag ? (SV_CALLING_SNIFFLES?.out?.vcf                 ?: null) : null
        sniffles_vcfgz                  = params.sniffles_flag ? (SV_CALLING_SNIFFLES?.out?.vcfgz               ?: null) : null
        sniffles_tbi                    = params.sniffles_flag ? (SV_CALLING_SNIFFLES?.out?.tbi                 ?: null) : null
        sniffles_csi                    = params.sniffles_flag ? (SV_CALLING_SNIFFLES?.out?.csi                 ?: null) : null

        sniffles_svync_vcf              = params.sniffles_flag ? (SV_CALLING_SNIFFLES?.out?.svync_vcf           ?: null) : null
        sniffles_svync_vcfgz            = params.sniffles_flag ? (SV_CALLING_SNIFFLES?.out?.svync_vcfgz         ?: null) : null
        sniffles_svync_tbi              = params.sniffles_flag ? (SV_CALLING_SNIFFLES?.out?.svync_tbi           ?: null) : null

        cutesv_vcfs                     = params.cutesv_flag ? (SV_CALLING_CUTESV?.out?.vcf                     ?: null) : null
        cutesv_vcfgz                    = params.cutesv_flag ? (SV_CALLING_CUTESV?.out?.vcfgz                   ?: null) : null
        cutesv_tbi                      = params.cutesv_flag ? (SV_CALLING_CUTESV?.out?.tbi                     ?: null) : null
        cutesv_csi                      = params.cutesv_flag ? (SV_CALLING_CUTESV?.out?.csi                     ?: null) : null

        cutesv_svync_vcf                = params.cutesv_flag ? (SV_CALLING_CUTESV?.out?.svync_vcf               ?: null) : null
        cutesv_svync_vcfgz              = params.cutesv_flag ? (SV_CALLING_CUTESV?.out?.svync_vcfgz             ?: null) : null
        cutesv_svync_tbi                = params.cutesv_flag ? (SV_CALLING_CUTESV?.out?.svync_tbi               ?: null) : null
}
