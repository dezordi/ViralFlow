process runIntraHostScript{
  publishDir "${params.outDir}/${sample_id}_results/"

  input:
     tuple val(sample_id), path(fa_bc), path(fa_algn)
     path(ref_gff)

  output:
     tuple val(sample_id), path("*.tsv"), path("*.fa")

  script:
     """
     python /intrahost_scriptv2.py \
            -in ${sample_id}.depth${params.depth}.fa.bc \
            -al ${sample_id}.depth${params.depth}.fa.algn \
            -dp ${params.minDpIntrahost} \
            -gf ${ref_gff}
     """
}
