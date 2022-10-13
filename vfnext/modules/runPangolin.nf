process runPangolin {
  publishDir "${params.outDir}/${sample_id}_results/",mode: "copy"
  input:
  tuple val(sample_id), path(intrahost_tsvs), path(algn_fasta), path(consensus_fa), path(ivar_txt), path(mut_tsv)
  //temporary solution, no need for ivar_txt and mut_tsv

  output:
   path("*.csv")
  shell:
  '''
  NUMLINES=$(wc -l < !{sample_id}.depth!{params.depth}.fa.bc.intrahost.short.tsv)

  if [ $NUMLINES -gt 1 ]; then
      cat !{sample_id}.depth!{params.depth}.fa !{sample_id}.depth!{params.depth}.fa.algn.minor.fa  > !{sample_id}.depth!{params.depth}.all.fa
      pangolin !{sample_id}.depth!{params.depth}.all.fa \
                -t !{params.pangolin_threads} --outfile !{sample_id}.all.fa.pango.out.csv
  fi
  if [ $NUMLINES -eq 1 ]; then
      pangolin !{sample_id}.depth!{params.depth}.fa \
                -t !{params.pangolin_threads} --outfile !{sample_id}.fa.pango.out.csv
  fi

  '''
}
