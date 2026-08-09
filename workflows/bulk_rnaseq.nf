/*
 * Bulk RNA-seq workflow.
 *
 * QC and trimming, three-way alignment/quantification (STAR, HISAT2, Salmon),
 * three-way differential expression (DESeq2, edgeR, limma-voom), KEGG/GO
 * enrichment, and cross-aligner / cross-method concordance reporting.
 */

include { FASTQC          } from '../modules/fastqc.nf'
include { TRIMGALORE      } from '../modules/trimgalore.nf'
include { STAR_ALIGN      } from '../modules/star.nf'
include { HISAT2_ALIGN    } from '../modules/hisat2.nf'
include { SALMON_QUANT    } from '../modules/salmon.nf'
include { FEATURECOUNTS as FEATURECOUNTS_STAR   } from '../modules/featurecounts.nf'
include { FEATURECOUNTS as FEATURECOUNTS_HISAT2 } from '../modules/featurecounts.nf'
include { MERGE_COUNTS as MERGE_STAR   } from '../modules/merge_counts.nf'
include { MERGE_COUNTS as MERGE_HISAT2 } from '../modules/merge_counts.nf'
include { TXIMPORT        } from '../modules/tximport.nf'
include { DESEQ2          } from '../modules/deseq2.nf'
include { EDGER           } from '../modules/edger.nf'
include { LIMMA_VOOM      } from '../modules/limma_voom.nf'
include { CLUSTERPROFILER } from '../modules/clusterprofiler.nf'
include { CONCORDANCE     } from '../modules/concordance.nf'

workflow BULK_RNASEQ {

    take:
    ch_samples
    refs_dir

    main:

    ch_versions = Channel.empty()

    // Quality control
    FASTQC(ch_samples)
    ch_versions = ch_versions.mix(FASTQC.out.versions)

    // Adapter and quality trimming
    TRIMGALORE(ch_samples)
    ch_versions = ch_versions.mix(TRIMGALORE.out.versions)

    ch_trimmed = TRIMGALORE.out.reads

    // ---- STAR alignment and gene-level quantification ---------------------
    STAR_ALIGN(ch_trimmed, "${refs_dir}/star_index")
    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions)

    // [] is the Nextflow idiom for an absent optional path.
    ch_star_bam = STAR_ALIGN.out.bam.map { sid, cond, bam ->
        tuple([id: sid, condition: cond, aligner: 'star',
               strandedness: params.strandedness, single_end: false], bam, [])
    }
    FEATURECOUNTS_STAR(ch_star_bam, "${refs_dir}/gencode.v44.annotation.gtf")
    ch_versions = ch_versions.mix(FEATURECOUNTS_STAR.out.versions)

    // Merge per-sample counts into one matrix. toSortedList keeps files and
    // IDs in the same order regardless of task completion order.
    ch_star_sorted = FEATURECOUNTS_STAR.out.counts.toSortedList { a, b -> a[0].id <=> b[0].id }
    MERGE_STAR(
        'star',
        ch_star_sorted.map { rows -> rows.collect { meta, f -> f } },
        ch_star_sorted.map { rows -> rows.collect { meta, f -> meta.id }.join(',') }
    )
    ch_versions = ch_versions.mix(MERGE_STAR.out.versions)

    // The DE modules take (meta, counts); meta.aligner drives their publishDir.
    ch_star_matrix = MERGE_STAR.out.counts.map { aligner, f -> tuple([aligner: aligner], f) }

    // ---- HISAT2 alignment and quantification ------------------------------
    HISAT2_ALIGN(ch_trimmed, "${refs_dir}/hisat2_index")
    ch_versions = ch_versions.mix(HISAT2_ALIGN.out.versions)

    ch_hisat2_bam = HISAT2_ALIGN.out.bam.map { sid, cond, bam, bai ->
        tuple([id: sid, condition: cond, aligner: 'hisat2',
               strandedness: params.strandedness, single_end: false], bam, bai)
    }

    FEATURECOUNTS_HISAT2(ch_hisat2_bam, "${refs_dir}/gencode.v44.annotation.gtf")
    ch_versions = ch_versions.mix(FEATURECOUNTS_HISAT2.out.versions)

    ch_hisat2_sorted = FEATURECOUNTS_HISAT2.out.counts.toSortedList { a, b -> a[0].id <=> b[0].id }
    MERGE_HISAT2(
        'hisat2',
        ch_hisat2_sorted.map { rows -> rows.collect { meta, f -> f } },
        ch_hisat2_sorted.map { rows -> rows.collect { meta, f -> meta.id }.join(',') }
    )
    ch_versions = ch_versions.mix(MERGE_HISAT2.out.versions)

    // ---- Salmon selective alignment with tximport aggregation -------------
    SALMON_QUANT(ch_trimmed, "${refs_dir}/salmon_index")
    ch_versions = ch_versions.mix(SALMON_QUANT.out.versions)

    TXIMPORT(
        SALMON_QUANT.out.quant.map { sid, cond, dir -> dir }.collect(),
        "${refs_dir}/tx2gene.tsv",
        file(params.input)
    )
    ch_versions = ch_versions.mix(TXIMPORT.out.versions)

    // ---- Differential expression on the merged STAR matrix ----------------

    ch_degs = Channel.empty()
    if (!params.skip_de) {

    DESEQ2(ch_star_matrix, file(params.input))
    EDGER(ch_star_matrix, file(params.input))
    LIMMA_VOOM(ch_star_matrix, file(params.input))
    ch_versions = ch_versions.mix(DESEQ2.out.versions, EDGER.out.versions, LIMMA_VOOM.out.versions)


    ch_all_degs = DESEQ2.out.results.map     { meta, f -> tuple(meta + [method: 'deseq2'], f) }
        .mix( EDGER.out.results.map          { meta, f -> tuple(meta + [method: 'edger'], f) } )
        .mix( LIMMA_VOOM.out.results.map     { meta, f -> tuple(meta + [method: 'limma_voom'], f) } )
    CLUSTERPROFILER(ch_all_degs)
    ch_versions = ch_versions.mix(CLUSTERPROFILER.out.versions)

    // ---- Cross-aligner and cross-method concordance (Figure 2A, 2B) -------
    CONCORDANCE(
        MERGE_STAR.out.counts.map   { aligner, f -> f },
        MERGE_HISAT2.out.counts.map { aligner, f -> f },
        TXIMPORT.out.counts,
        DESEQ2.out.results.map      { meta, f -> f },
        EDGER.out.results.map       { meta, f -> f },
        LIMMA_VOOM.out.results.map  { meta, f -> f }
    )
    ch_versions = ch_versions.mix(CONCORDANCE.out.versions)

    ch_degs = DESEQ2.out.results

    }   // end if (!params.skip_de)

    emit:
    counts   = MERGE_STAR.out.counts.map { aligner, f -> f }
    degs     = ch_degs
    versions = ch_versions
}
