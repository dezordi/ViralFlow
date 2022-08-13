process runNextClade {
  publishDir "${params.outDir}/${sample_id}_results/"

  input:
  tuple val(sample_id), path(intrahost_tsvs), path(algn_fasta), path(consensus_fa), path(ivar_txt), path(mut_tsv)
  path(ref_fa)

  // temporary solution, no need for ivar_txt and mut_tsv
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
      nextclade run --jobs !{params.threads} \
                --input-root-seq=!{ref_fa} \
                --input-dataset=!{nxt_dataset} \
                --output-csv=!{sample_id}.depth!{params.depth}.all.fa.nextclade.csv \
                --output-all=./ \
                !{sample_id}.depth!{params.depth}.all.fa
  fi

  if [ $NUMLINES -eq 1 ]; then
      nextclade run --jobs !{params.threads} \
             --input-root-seq=!{ref_fa} \
             --input-dataset=!{nxt_dataset} \
             --output-csv=!{sample_id}.depth!{params.depth}.fa.nextclade.csv \
             --output-all=./ \
             !{sample_id}.depth!{params.depth}.fa
  fi

  '''
}
