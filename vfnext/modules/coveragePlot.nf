process coveragePlot {
    publishDir "${params.outDir}/${sample_id}_results/", mode: "copy"

    input:

      tuple val(sample_id), path(bam_files), path(bai_files)
      val(genome_code)
      path(refIndexFiles)
    output:
        path("*coveragePlot*")

    script:
    
    bam = bam_files[0].toString() // Certifique-se de que o índice [0] seja apropriado para o seu caso
    genomecode = genome_code.toString()
    depth = params.depth

    """
    bamdash -b ${bam} -r ${genomecode} -c ${depth}
    mv ${genomecode}_plot.html ${sample_id}_coveragePlot.html
    bamdash -b ${bam} -r ${genomecode} -c ${depth} -e svg
    mv ${genomecode}_plot.svg ${sample_id}_coveragePlot.svg  
    """
}
