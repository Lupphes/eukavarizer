/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES FOR STRUCTURAL VARIANT MERGING AND SYNCHRONIZATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SURVIVOR_MERGE        } from '../../../modules/nf-core/survivor/merge/main'
include { SURVIVOR_STATS        } from '../../../modules/nf-core/survivor/stats/main'
include { BCFTOOLS_MERGE        } from '../../../modules/nf-core/bcftools/merge/main'
include { VCF_REPGEN            } from '../../../modules/local/vcf_repgen/main'

workflow SV_UNIFICATION {

    // take:



    // main:

}
