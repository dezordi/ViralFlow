process runNextClade {
  publishDir "${params.outDir}/${sample_id}_results/"

  input:
  tuple val(sample_id), path(intrahost_tsvs), path(algn_fasta), path(consensus_fa)

  output:
  tuple path("*.csv"), path("*.fasta")

  shell:
  nxt_dataset = "${workflow.projectDir}/containers/nextclade_dataset/sars-cov-2/"
  '''
  # check if intravariants are present:
  #     if yes, compile all on a single file and use it as input
  #     if no, use the consensus sequence as input

  NUMLINES=$(wc -l < !{sample_id}.depth!{params.depth}.fa.bc.intrahost.short.tsv)

  if [ $NUMLINES -gt 1 ]; then
      cat !{sample_id}.depth!{params.depth}.fa !{sample_id}.depth!{params.depth}.fa.algn.minor.fa  > !{sample_id}.depth!{params.depth}.all.fa
      nextclade -i !{sample_id}.depth!{params.depth}.all.fa --jobs !{params.threads} \
                --input-root-seq=!{params.referenceGenome} \
                --input-dataset=!{nxt_dataset} \
                --output-csv=!{sample_id}.depth!{params.depth}.all.fa.nextclade.csv \
                --output-dir=./
  fi

  if [ $NUMLINES -eq 1 ]; then
      nextclade -i !{sample_id}.depth!{params.depth}.fa --jobs !{params.threads} \
             --input-root-seq=!{params.referenceGenome} \
             --input-dataset=!{nxt_dataset} \
             --output-csv=!{sample_id}.depth!{params.depth}.fa.nextclade.csv \
             --output-dir=./
  fi

  '''

}
