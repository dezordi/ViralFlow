process coveragePlot {
    publishDir "${params.outDir}/${sample_id}_results/", mode: "copy"

    input:

      tuple val(sample_id), path(bam_files), path(bai_files)
      val(genome_code)
    
    output:
        path("*coveragePlot*"), optional: true
        path("coveragePlot_result.txt"), optional: true, emit: result, hidden: true
    script:
    
    bam = bam_files[0].toString()
    genomecode = genome_code.toString()
    depth = params.depth
    html = "${sample_id}_coveragePlot.html"
    png = "${sample_id}_coveragePlot.png"
    svg = "${sample_id}_coveragePlot.svg"
    
    """
    #!/usr/bin/env python
    # ----- import libraries -------------------------------------------------
    import pysam
    import subprocess
    
    n_mapped_reads = pysam.AlignmentFile('${bam}','rb').count()
    if n_mapped_reads > 0:
      subprocess.run([f"bamdash -b ${bam} -r ${genomecode} -c ${depth} -e svg"], shell=True)
      subprocess.run([f"mv ${genomecode}_plot.svg ${svg}"], shell=True)

      subprocess.run([f"bamdash -b ${bam} -r ${genomecode} -c ${depth} -e png"], shell=True)
      subprocess.run([f"mv ${genomecode}_plot.png ${png}"], shell=True)

      subprocess.run([f"bamdash -b ${bam} -r ${genomecode} -c ${depth}"], shell=True)
      subprocess.run([f"mv ${genomecode}_plot.html ${html}"], shell=True)
   
    else:
      result = "No mapped reads were found in the sorted BAM file for sample ${sample_id}. The coverage plot will not be generated for it."
      with open('coveragePlot_result.txt', 'w') as f:
        f.write(result)
      """
}

process snpPlot {
    publishDir "${params.outDir}/${sample_id}_results/", mode: "copy"

    input:

      tuple val(sample_id), path("${sample_id}.depth${params.depth}.fa.algn")

      
    output:
        path("*snpPlot*")

    script:
    
    plot_name = "${sample_id}_snpPlot"
    

    """

    snipit ${sample_id}.depth${params.depth}.fa.algn -o ${plot_name} --solid-background 
    snipit ${sample_id}.depth${params.depth}.fa.algn -o ${plot_name} -f svg --solid-background

    """
}
