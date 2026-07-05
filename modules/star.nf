/*
 * STAR v2.7.11a
 *
 * Two-pass alignment against GRCh38.p14 / GENCODE v44.
 * STAR is the primary alignment-based aligner in KazRNA-Pipe; its output
 * feeds the primary count matrix used for downstream DE and deconvolution.
 */

process STAR_ALIGN {
    tag "${sample_id}"
    label 'process_high'
    publishDir "${params.outdir}/bulk/star", mode: params.publish_dir_mode

    container 'quay.io/biocontainers/star:2.7.11a--h0033a41_0'

    input:
    tuple val(sample_id), val(condition), path(fastq_1), path(fastq_2)
    path  star_index

    output:
    tuple val(sample_id), val(condition),
          path("${sample_id}.Aligned.sortedByCoord.out.bam"),
          emit: bam
    path "${sample_id}.Log.final.out",        emit: log
    path "${sample_id}.SJ.out.tab",           emit: sj
    path "${sample_id}.timing.tsv",           emit: timing
    path 'versions.yml',                      emit: versions

    script:
    """
    /usr/bin/time -f '%e\\t%M' -o ${sample_id}.timing.raw \\
    STAR \\
        --runThreadN ${task.cpus} \\
        --genomeDir ${star_index} \\
        --readFilesIn ${fastq_1} ${fastq_2} \\
        --readFilesCommand zcat \\
        --twopassMode Basic \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMstrandField intronMotif \\
        --quantMode GeneCounts \\
        --outFileNamePrefix ${sample_id}. \\
        --limitBAMsortRAM 30000000000

    # Convert /usr/bin/time output (seconds <TAB> kbytes_max_rss) to a tidy TSV
    awk -v s=${sample_id} 'BEGIN{print "sample\\tstep\\twall_seconds\\tmax_rss_kb"}
        {print s "\\tstar\\t" \$1 "\\t" \$2}' ${sample_id}.timing.raw > ${sample_id}.timing.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$( STAR --version )
    END_VERSIONS
    """
}
