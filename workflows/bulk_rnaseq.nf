/*
 * Bulk RNA-seq workflow.
 *
 * QC and trimming, three-way alignment/quantification (STAR, HISAT2, Salmon),
 * three-way differential expression (DESeq2, edgeR, limma-voom), KEGG/GO
 * enrichment, and cross-aligner / cross-method concordance reporting (Figure 2).
 *
 * A single `meta` map ([id, condition]) is carried through every channel so the
 * module interfaces remain consistent from alignment to reporting.
 */

include { FASTQC          } from '../modules/fastqc.nf'
include { TRIMGALORE      } from '../modules/trimgalore.nf'
include { STAR_ALIGN      } from '../modules/star.nf'
include { HISAT2_ALIGN    } from '../modules/hisat2.nf'
include { SALMON_QUANT    } from '../modules/salmon.nf'
include { FEATURECOUNTS   } from '../modules/featurecounts.nf'
include { TXIMPORT        } from '../modules/tximport.nf'
include { DESEQ2          } from '../modules/deseq2.nf'
include { EDGER           } from '../modules/edger.nf'
include { LIMMA_VOOM      } from '../modules/limma_voom.nf'
include { CLUSTERPROFILER } from '../modules/clusterprofiler.nf'
include { CONCORDANCE     } from '../modules/concordance.nf'

workflow BULK_RNASEQ {

    take:
    ch_samples      // (sample_id, condition, fastq_1, fastq_2)
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

    // STAR alignment and gene-level quantification
    STAR_ALIGN(ch_trimmed, "${refs_dir}/star_index")
    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions)

    ch_star_bam = STAR_ALIGN.out.bam.map { sid, cond, bam ->
        tuple([id: sid, condition: cond], bam, file('NO_BAI'))
    }
    FEATURECOUNTS(ch_star_bam, "${refs_dir}/gencode.v44.annotation.gtf")
    ch_versions = ch_versions.mix(FEATURECOUNTS.out.versions)

    ch_star_counts    = FEATURECOUNTS.out.counts
    ch_star_collected = ch_star_counts.map { meta, f -> f }.collect()

    // HISAT2 alignment and quantification
    HISAT2_ALIGN(ch_trimmed, "${refs_dir}/hisat2_index")
    ch_versions = ch_versions.mix(HISAT2_ALIGN.out.versions)
    ch_hisat2_collected = HISAT2_ALIGN.out.counts.collect()

    // Salmon selective alignment with tximport gene-level aggregation
    SALMON_QUANT(ch_trimmed, "${refs_dir}/salmon_index")
    ch_versions = ch_versions.mix(SALMON_QUANT.out.versions)
    TXIMPORT(
        SALMON_QUANT.out.quant.map { sid, cond, q -> q }.collect(),
        "${refs_dir}/tx2gene.tsv",
        params.input
    )
    ch_salmon_collected = TXIMPORT.out.counts
    ch_versions = ch_versions.mix(TXIMPORT.out.versions)

    // Gene-level counts from all three aligners for concordance analysis
    ch_counts_all = ch_star_collected.map   { tuple('star',   it) }
        .concat( ch_hisat2_collected.map { tuple('hisat2', it) } )
        .concat( ch_salmon_collected.map { tuple('salmon', it) } )

    // Differential expression with three methods on the primary STAR counts
    DESEQ2(ch_star_counts, params.input)
    EDGER(ch_star_counts, params.input)
    LIMMA_VOOM(ch_star_counts, params.input)
    ch_versions = ch_versions.mix(DESEQ2.out.versions, EDGER.out.versions, LIMMA_VOOM.out.versions)

    // KEGG and GO enrichment on each differential expression result
    ch_all_degs = DESEQ2.out.results
        .mix(EDGER.out.results, LIMMA_VOOM.out.results)
    CLUSTERPROFILER(ch_all_degs)
    ch_versions = ch_versions.mix(CLUSTERPROFILER.out.versions)

    // Cross-aligner and cross-method concordance (Figure 2A, 2B)
    ch_de_all = DESEQ2.out.results.map   { meta, f -> f }
        .mix( EDGER.out.results.map      { meta, f -> f } )
        .mix( LIMMA_VOOM.out.results.map { meta, f -> f } )
        .collect()
    CONCORDANCE(
        ch_counts_all.map { label, f -> f }.collect(),
        ch_de_all,
        params.input
    )
    ch_versions = ch_versions.mix(CONCORDANCE.out.versions)

    emit:
    counts   = ch_star_counts
    degs     = DESEQ2.out.results
    versions = ch_versions
}
