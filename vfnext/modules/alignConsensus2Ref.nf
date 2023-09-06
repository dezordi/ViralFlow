process alignConsensus2Ref {
    publishDir "${params.outDir}/${sample_id}_results/", mode: "copy"
    label "multithread"

    input:
    tuple val(sample_id), path(consensus_fa), path(ivar_txt), path(mut_tsv)
    path(ref_fa)
    //temporary solution, we only need the consensus_fa. handle channels better
    output:
    tuple val(sample_id), path("${sample_id}.depth${params.depth}.fa.algn")
    //path("*.fa.algn")
    script:
    """
    mafft --keeplength --add ${sample_id}.depth${params.depth}.fa \
                       --thread ${params.mafft_threads} \
                       ${ref_fa} \
    > ${sample_id}.depth${params.depth}.fa.algn
    """
}
