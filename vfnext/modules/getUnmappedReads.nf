process getUnmappedReads {
    //publishDir "${params.outDir}/${sample_id}_results/"
    label "singlethread"
    input:
        tuple val(sample_id), path(bam_files)

    output:
        tuple val(sample_id), path("unmapped.bam")

    script:
    """
    samtools view -b -f 4 ${sample_id}.sorted.bam > unmapped.bam
    samtools sort -n -o unmapped.bam unmapped.bam
#    samtools view -hbf 64 unmapped.bam > unmapped.R1.bam
#    samtools view -hbf 128 unmapped.bam > unmapped.R2.bam
    """
}
