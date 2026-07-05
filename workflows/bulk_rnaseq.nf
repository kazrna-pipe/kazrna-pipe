/*
 * Bulk RNA-seq workflow.
 *
 * QC → trim → 3-way alignment/quant → 3-way DE → KEGG enrichment → concordance.
 * Three aligners and three DE methods run in parallel so that the final stage
 * can report cross-method concordance (manuscript Figure 2).
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

    // 1. QC
    FASTQC(ch_samples)
    ch_versions = ch_versions.mix(FASTQC.out.versions)

    // 2. Adapter / quality trim
    TRIMGALORE(ch_samples)
    ch_versions = ch_versions.mix(TRIMGALORE.out.versions)

    ch_trimmed = TRIMGALORE.out.reads

    // 3a. STAR (alignment-based)
    STAR_ALIGN(ch_trimmed, "${refs_dir}/star_index")
    FEATURECOUNTS(STAR_ALIGN.out.bam, "${refs_dir}/gencode.v44.annotation.gtf")
    ch_star_counts = FEATURECOUNTS.out.counts.collect()
    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions, FEATURECOUNTS.out.versions)

    // 3b. HISAT2 (alignment-based)
    HISAT2_ALIGN(ch_trimmed, "${refs_dir}/hisat2_index")
    ch_hisat2_counts = HISAT2_ALIGN.out.counts.collect()
    ch_versions = ch_versions.mix(HISAT2_ALIGN.out.versions)

    // 3c. Salmon (selective alignment)
    SALMON_QUANT(ch_trimmed, "${refs_dir}/salmon_index")
    TXIMPORT(SALMON_QUANT.out.quant.collect(), "${refs_dir}/tx2gene.tsv")
    ch_salmon_counts = TXIMPORT.out.counts
    ch_versions = ch_versions.mix(SALMON_QUANT.out.versions, TXIMPORT.out.versions)

    // Merge the three count matrices into one channel keyed by aligner
    ch_counts_all = Channel.empty()
        .mix( ch_star_counts.map   { tuple('star',   it) } )
        .mix( ch_hisat2_counts.map { tuple('hisat2', it) } )
        .mix( ch_salmon_counts.map { tuple('salmon', it) } )

    // 4. Differential expression - run on STAR counts for the primary analysis
    DESEQ2(ch_star_counts, ch_samples.map { it[0..1] }.collect())
    EDGER(ch_star_counts, ch_samples.map { it[0..1] }.collect())
    LIMMA_VOOM(ch_star_counts, ch_samples.map { it[0..1] }.collect())
    ch_versions = ch_versions.mix(DESEQ2.out.versions, EDGER.out.versions, LIMMA_VOOM.out.versions)

    // 5. Pathway enrichment on the core (3-way intersection) DEG set
    CLUSTERPROFILER(
        DESEQ2.out.degs,
        EDGER.out.degs,
        LIMMA_VOOM.out.degs
    )
    ch_versions = ch_versions.mix(CLUSTERPROFILER.out.versions)

    // 6. Cross-aligner and cross-method concordance reports (Figure 2A, 2B)
    CONCORDANCE(
        ch_counts_all.collect(),
        DESEQ2.out.degs,
        EDGER.out.degs,
        LIMMA_VOOM.out.degs
    )
    ch_versions = ch_versions.mix(CONCORDANCE.out.versions)

    emit:
    counts   = ch_star_counts                 // primary count matrix for downstream deconv
    degs     = DESEQ2.out.degs                // primary DEG table
    versions = ch_versions
}
