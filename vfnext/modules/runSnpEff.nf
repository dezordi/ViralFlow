process runSnpEff{
  errorStrategy 'ignore'
  publishDir "${params.outDir}/${sample_id}_results/"

  input:
    tuple val(sample_id), path(bam_files)
    val(genome_code)
    path(refGenomeFasta)
    path(refIndexFiles)
    path(entry_found_file) // here to assure runSnpEff will not start before DB was checked
  output:
  tuple val(sample_id), path("*.vcf"), path("snpEff_summary.html"), path("snpEff_genes.txt")

  script:
  ref_fa = "${refGenomeFasta}"
  sorted_bam = "${bam_files[0].getSimpleName()}.sorted.bam"
  """
  freebayes -p 1 --reference-quality 30 \
            -f ${ref_fa} ${sorted_bam} > ${sample_id}.vcf
  snpEff download -v ${genome_code}
  snpEff ann -Xmx4g \
            ${genome_code} ${sample_id}.vcf > ${sample_id}.ann.vcf
  """
}
