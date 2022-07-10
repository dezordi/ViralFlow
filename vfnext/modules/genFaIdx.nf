process genFaIdx {
    /*
    * Indexes reference fasta file using bwa.
    */
    publishDir "${params.outDir}/"

    input:
        path(reference_fasta)

    output:
        path("${reference_fasta}*")

    script:
        """
        samtools faidx ${reference_fasta}
        """
}
