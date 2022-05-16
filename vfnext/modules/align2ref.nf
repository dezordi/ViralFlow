process align2ref{
  input:
    tuple val(sample_id), path(reads), path(fasta_amb), path(fasta_ann), path(fasta_bwt), path(fasta_pac), path(fasta_sa)
//    path("*.fasta*")

  output:
    tuple val(sample_id), path("*.bam")


  script:
    """
    # run aligner -------------------------------------------------------------
    ln -s ${params.referenceGenome} ./${fasta_amb.getSimpleName()}.fasta
    bwa mem -t ${params.threads} ./${fasta_amb.getSimpleName()}.fasta ${reads[0]} ${reads[1]} \
            -o ${sample_id}.bam

    # Sort alignments by leftmost coordinates ---------------------------------
    samtools sort -o ${sample_id}.sorted.bam ${sample_id}.bam

    # Index a coordinate-sorted BGZIP-compressed SAM, BAM or CRAM file for fast
    # random access.
    samtools index ${sample_id}.sorted.bam
    """
}
