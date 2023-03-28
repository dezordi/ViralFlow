
process align2ref{
  publishDir "${params.outDir}/${sample_id}_results/", mode: "copy"
  
  input:
    tuple val(sample_id), path(reads), path(fasta_amb), path(fasta_ann), path(fasta_bwt), path(fasta_pac), path(fasta_sa)
    path(ref_fa)

  output:
    tuple val(sample_id), path("*.sorted.bam"), path("*.bai"), emit: regular_output
    path("${sample_id}.trimmed_reads.txt"), emit:trimmed_reads, optional: true
  script:
    trim_bam = "${sample_id}.trimmed"
    bed = "${params.primersBED}"
    if (params.primersBED==null)
    """
    # run aligner -------------------------------------------------------------
    #ln -s ${ref_fa} ./${fasta_amb.getSimpleName()}.fasta
    bwa mem ./${ref_fa} ${reads[0]} ${reads[1]} \
            -o ${sample_id}.bam -t ${params.bwa_threads}

    # Sort alignments by leftmost coordinates ---------------------------------
    samtools sort -o ${sample_id}.sorted.bam ${sample_id}.bam

    # Index a coordinate-sorted BGZIP-compressed SAM, BAM or CRAM file for fast
    # random access.
    samtools index ${sample_id}.sorted.bam    
    """

    else if (!(params.primersBED==null))
    """
    # run aligner -------------------------------------------------------------
    #ln -s ${ref_fa} ./${fasta_amb.getSimpleName()}.fasta
    bwa mem ./${ref_fa} ${reads[0]} ${reads[1]} \
            -o ${sample_id}.bam -t ${params.bwa_threads}

    # Sort alignments by leftmost coordinates ---------------------------------
    samtools sort -o ${sample_id}.sorted.bam ${sample_id}.bam

    # Index a coordinate-sorted BGZIP-compressed SAM, BAM or CRAM file for fast
    # random access.
    samtools index ${sample_id}.sorted.bam    
    
    #TRIM PRIMERS OF SORTED BAM

    samtools ampliconclip --both-ends --hard-clip --filter-len ${params.minLen} -b ${bed} ${sample_id}.sorted.bam -f ${sample_id}.trimmed_reads.txt > ${trim_bam}
    # Sort alignments by leftmost coordinates ---------------------------------
    samtools sort -o ${trim_bam}.sorted.bam ${trim_bam}
    samtools index ${trim_bam}.sorted.bam
    mv ${sample_id}.sorted.bam ${sample_id}.raw.sorted.bam
    mv ${sample_id}.sorted.bam.bai ${sample_id}.raw.sorted.bam.bai
    mv ${trim_bam}.sorted.bam ${sample_id}.sorted.bam
    mv ${trim_bam}.sorted.bam.bai ${sample_id}.sorted.bam.bai
    """
}
