
process runFastp{
  publishDir "${params.outDir}/${sample_id}_results/"
  label "multithread"

  input:
    tuple val(sample_id), path(reads)

  output:
    tuple val(prfx), path('*.R{1,2}.fq.gz'), path("${prfx}.fastp.html")

  script:
      prfx = reads[0].getSimpleName()
      if (params.adaptersFile==null)
      """
      fastp -i ${reads[0]} -I ${reads[1]} \
            --thread ${params.fastp_threads} \
            -o ${reads[0].getSimpleName()}.R1.fq.gz \
            -O ${reads[0].getSimpleName()}.R2.fq.gz \
            -h ${reads[0].getSimpleName()}.fastp.html \
            -j ${reads[0].getSimpleName()}.fastp.json \
            -l ${params.minLen} -f ${params.trimLen} -t ${params.trimLen} \
            -F ${params.trimLen} -T ${params.trimLen} \
            --cut_front --cut_tail --qualified_quality_phred 20
      """

      else if (!(params.adaptersFile==null))
      """
      fastp -i ${reads[0]} -I ${reads[1]} \
            --adapter_fasta ${params.adaptersFile} \
            --thread ${params.fastp_threads} \
            -o ${reads[0].getSimpleName()}.R1.fq.gz -O ${reads[0].getSimpleName()}.R2.fq.gz -h ${reads[0].getSimpleName()}.fastp.html \
            -j ${reads[0].getSimpleName()}.fastp.json \
            -l ${params.minLen} -f ${params.trimLen} -t ${params.trimLen} \
            -F ${params.trimLen} -T ${params.trimLen} \
            --cut_front --cut_tail --qualified_quality_phred 20
      """
}
