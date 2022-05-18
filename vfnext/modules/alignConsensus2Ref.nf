process alignConsensus2Ref {

    input:
    tuple val(sample_id), path(consensus_fa)

    //output:
    //path("${sample_id}.depth${depth}.fa.algn")

    script:
    """
    mafft --keeplength --add ${sample_id}.depth${params.depth}.fa \
                       --thread ${params.threads} \
                       ${params.referenceGenome} \
    > ${sample_id}.depth${params.depth}.fa.algn
    """
}
