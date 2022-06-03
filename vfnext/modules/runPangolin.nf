process runPangolin {
  input:
  tuple val(sample_id), path(intrahost_tsvs), path(algn_fasta), path(consensus_fa)

  shell:
  '''
  NUMLINES=$(wc -l < !{sample_id}.depth!{params.depth}.fa.bc.intrahost.short.tsv)

  if [ $NUMLINES -gt 1 ]; then
      cat !{sample_id}.depth!{params.depth}.fa !{sample_id}.depth!{params.depth}.fa.algn.minor.fa  > !{sample_id}.depth!{params.depth}.all.fa
      pangolin !{sample_id}.depth!{params.depth}.all.fa \
                -t !{params.threads} --outfile !{sample_id}.all.fa.pango.out_csv
  fi
  if [ $NUMLINES -eq 1 ]; then
      pangolin !{sample_id}.depth!{params.depth}.fa \
                -t !{params.threads} --outfile !{sample_id}.fa.pango.csv
  fi

  '''
}
