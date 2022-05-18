process runReadCounts{
  input:
  tuple val(sample_id), path(bams)

  output:
  path("${sample_id}.depth${d}.fa.bc")

  script:
  sorted_bam = "${bams[0].getSimpleName()}.sorted.bam"
  ref_gnm = "${params.referenceGenome}"
  d = "${params.depth}"
  """
  # RUN READ COUNT
  bam-readcount -d 50000 -q 30 -w 0 -f ${ref_gnm} ${sorted_bam} > ${sample_id}.depth${d}.fa.bc
  """
}
