process runIntraHostScript{

  input:
     tuple val(sample_id), path(fa_bc), path(fa_algn)

  script:
     """
     python /intrahost_script.py \
            -in ${sample_id}.depth${params.depth}.fa.bc \
            -al ${sample_id}.depth${params.depth}.fa.algn \
            -dp ${params.minDpIntrahost}
     """
}
