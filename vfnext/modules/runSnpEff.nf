process runSnpEff{
  errorStrategy 'ignore'
  publishDir "${params.outDir}/${sample_id}_results/"

  input:
  tuple val(sample_id), path(bam_files)
  val(genome_code)
  path(refGenomeFasta)
  path(refIndexFiles)
  output:
  tuple val(sample_id), path("*.vcf")

  script:
  ref_fa = "${refGenomeFasta}"
  sorted_bam = "${bam_files[0].getSimpleName()}.sorted.bam"
  """
  freebayes -p 1 --reference-quality 30 \
            -f ${ref_fa} ${sorted_bam} > ${sample_id}.vcf
  snpEff download -v ${genome_code}
  snpEff ann -Xmx4g -noStats \
            ${genome_code} ${sample_id}.vcf > ${sample_id}.ann.vcf
  """
}
