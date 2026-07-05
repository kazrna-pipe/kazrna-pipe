nextflow.enable.dsl = 2

include { SEURAT_WORKFLOW } from '../modules/seurat.nf'

workflow SCRNASEQ_CPU {
    take:
    sc_samples      

    main:
    SEURAT_WORKFLOW(sc_samples)

    emit:
    reference = SEURAT_WORKFLOW.out.reference 
    versions  = SEURAT_WORKFLOW.out.versions
}