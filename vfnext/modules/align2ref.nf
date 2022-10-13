
process align2ref{
  publishDir "${params.outDir}/${sample_id}_results/", mode: "copy", pattern: "*.{sorted.bam, bai}"
  
  input:
    tuple val(sample_id), path(reads), path(fasta_amb), path(fasta_ann), path(fasta_bwt), path(fasta_pac), path(fasta_sa)
    path(ref_fa)

  output:
    tuple val(sample_id), path("*.sorted.bam"), path("*.bai")

  script:
    """
    # run aligner -------------------------------------------------------------
    #ln -s ${ref_fa} ./${fasta_amb.getSimpleName()}.fasta
    bwa mem -t ${params.threads} ./${ref_fa} ${reads[0]} ${reads[1]} \
            -o ${sample_id}.bam -t ${params.bwa_threads}

    # Sort alignments by leftmost coordinates ---------------------------------
    samtools sort -o ${sample_id}.sorted.bam ${sample_id}.bam

    # Index a coordinate-sorted BGZIP-compressed SAM, BAM or CRAM file for fast
    # random access.
    samtools index ${sample_id}.sorted.bam
    """
}
