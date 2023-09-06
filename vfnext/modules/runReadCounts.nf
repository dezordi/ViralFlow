process runReadCounts{
  publishDir "${params.outDir}/${sample_id}_results/", mode: "copy"

  input:
  tuple val(sample_id), path(bams)
  path(ref_fa)

  output:
  tuple val(sample_id), path("${sample_id}.depth${d}.fa.bc")

  script:
  sorted_bam = "${bams[0].getSimpleName()}.sorted.bam"
  d = "${params.depth}"
  """
  # RUN READ COUNT
  bam-readcount -d 50000 -q 30 -w 0 -f ${ref_fa} ${sorted_bam} > ${sample_id}.depth${d}.fa.bc
  """
}
