process compileOutputs{
  publishDir "${params.outDir}/"
  input:
    val(go)

  output:
    path("*")
  script:
    """
    compileOutput.py -dD ${params.outDir} \
                            -oD ./ \
                            --depth ${params.depth}
    """
}
