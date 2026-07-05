nextflow.enable.dsl = 2

include { RAPIDS_SINGLECELL } from '../modules/rapids_singlecell.nf'

workflow SCRNASEQ_GPU {
    take:
    sc_samples        

    main:
    RAPIDS_SINGLECELL(sc_samples)

    emit:
    versions = RAPIDS_SINGLECELL.out.versions
}