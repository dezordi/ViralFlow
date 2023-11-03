process coveragePlot {
    publishDir "${params.outDir}/${sample_id}_results/", mode: "copy"

    input:

      tuple val(sample_id), path(bam_files), path(bai_files)
      val(genome_code)
    output:
        path("*coveragePlot*")

    script:
    
    bam = bam_files[0].toString()
    genomecode = genome_code.toString()
    depth = params.depth
    html = "${sample_id}_coveragePlot.html"
    png = "${sample_id}_coveragePlot.png"
    svg = "${sample_id}_coveragePlot.svg"

    """
    bamdash -b ${bam} -r ${genomecode} -c ${depth} -e svg
    mv ${genomecode}_plot.svg ${svg}

    bamdash -b ${bam} -r ${genomecode} -c ${depth} -e png
    mv ${genomecode}_plot.png ${png}

    bamdash -b ${bam} -r ${genomecode} -c ${depth}
    mv ${genomecode}_plot.html ${html}
    """
}
