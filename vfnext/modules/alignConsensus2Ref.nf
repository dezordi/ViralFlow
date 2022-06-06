process alignConsensus2Ref {
   publishDir "${params.outDir}/${sample_id}_results/"
    input:
    tuple val(sample_id), path(consensus_fa)

    output:
    tuple val(sample_id), path("${sample_id}.depth${params.depth}.fa.algn")
    //path("*.fa.algn")
    script:
    """
    mafft --keeplength --add ${sample_id}.depth${params.depth}.fa \
                       --thread ${params.threads} \
                       ${params.referenceGenome} \
    > ${sample_id}.depth${params.depth}.fa.algn
    """
}
