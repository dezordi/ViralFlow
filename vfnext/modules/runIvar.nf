process runIvar{

  input:
    tuple val(sample_id), path(bams)

  output:
    tuple val(sample_id), path("*.depth*.fa")

  script:
    sorted_bam = "${bams[0].getSimpleName()}.sorted.bam"
    ref_gnm = "${params.referenceGenome}"
    d = "${params.depth}"
    """
    # IVAR STEP 1 ----------------------------------------------------------------
    samtools mpileup -aa -d 50000 --reference ${ref_gnm} -a -B ${sorted_bam} | \
       ivar variants -p ${sample_id} -q 30 -t 0.05

    # IVAR STEP 2 ----------------------------------------------------------------
    samtools mpileup -aa -d 50000 --reference ${ref_gnm} -a -B ${sorted_bam} | \
       ivar consensus -p ${sample_id} -q 30 -t 0 -m ${d} -n N

    # IVAR STEP 3 ----------------------------------------------------------------
    samtools mpileup -aa -d 50000 --reference ${ref_gnm} -a -B -B ${sorted_bam} | \
       ivar consensus -p ${sample_id}.ivar060 -q 30 -t 0.60 -n N -m ${params.depth}
    # EDIT FILE NAMES
    mv ${sample_id}.fa ${sample_id}.depth${d}.fa
    mv ${sample_id}.ivar060.fa ${sample_id}.depth${d}.amb.fa
    sed -i -e 's/>.*/>${sample_id}/g' ${sample_id}.depth${d}.fa
    sed -i -e 's/>.*/>${sample_id}/g' ${sample_id}.depth${d}.amb.fa
    """
}
