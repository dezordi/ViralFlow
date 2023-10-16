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
    """
}
