process compileOutputs{
  publishDir "${params.outDir}/COMPILED_OUTPUT/"
  label "singlethread"
  
  input:
    val(go)
    val(virus_tag)

  output:
    path("*")
  script:
    """
    compileOutput.py -dD ${params.outDir} \
                            -oD ./ \
                            --depth ${params.depth} \
                            -virus_tag ${virus_tag}
    """
}
