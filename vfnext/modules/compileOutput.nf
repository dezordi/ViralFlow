process compileOutputs{
  publishDir "${params.outDir}/COMPILED_OUTPUT/"
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
