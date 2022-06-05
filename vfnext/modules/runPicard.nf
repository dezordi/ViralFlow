process runPicard {
  input:
    tuple val(sample_id), path(bams)

  output:
    tuple val(sample_id), path("wgs")

  script:
  sorted_bam = "${sample_id}.sorted.bam"
  java_cmd = "java -jar /app/picard.jar"
  """
  ${java_cmd} CollectWgsMetrics -I ${sorted_bam} \
                                -R ${params.referenceGenome}\
                                -O wgs -CAP 99999

  ${java_cmd} CollectMultipleMetrics -I ${sorted_bam} \
                                     -R ${params.referenceGenome}\
                                     -O metrics
  """
}
