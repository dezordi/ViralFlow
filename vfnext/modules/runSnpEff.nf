process runSnpEff{
  errorStrategy 'ignore'
  publishDir "${params.outDir}/${sample_id}_results/", mode: "copy"

  input:
    tuple val(sample_id), path(bam_files)
    val(genome_code)
    path(refGenomeFasta)
    path(refIndexFiles)
    path(entry_found_file) // here to assure runSnpEff will not start before DB was checked
  
  output:
  tuple val(sample_id), path("*.vcf"), path("${sample_id}_snpEff_summary.html"), path("snpEff_genes.txt")

  script:
  ref_fa = "${refGenomeFasta}"
  sorted_bam = "${bam_files[0].getSimpleName()}.sorted.bam"
  if (params.virus == "custom")
    """
    freebayes -p 1 --reference-quality 30 \
            -f ${ref_fa} ${sorted_bam} > ${sample_id}.vcf
    snpEff ann -Xmx4g \
            -v ${genome_code} ${sample_id}.vcf > ${sample_id}.ann.vcf
    # add sample id to htmls
    mv snpEff_summary.html ${sample_id}_snpEff_summary.html
    """
  else
    """
    freebayes -p 1 --reference-quality 30 \
            -f ${ref_fa} ${sorted_bam} > ${sample_id}.vcf
    snpEff download -v ${genome_code}
    snpEff ann -Xmx4g \
            ${genome_code} ${sample_id}.vcf > ${sample_id}.ann.vcf
    # add sample id to htmls
    mv snpEff_summary.html ${sample_id}_snpEff_summary.html
    """
}
