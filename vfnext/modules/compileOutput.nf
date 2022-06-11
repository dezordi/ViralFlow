process compileOutputs{
  input:
    val(go)
  script:
    """
    python /compileOutput.py -dD ${params.outDir} -oD ${params.outDir} \
                             --depth ${params.depth}
    """
}
