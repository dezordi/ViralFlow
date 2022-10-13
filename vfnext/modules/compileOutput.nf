process compileOutputs{
  publishDir "${params.outDir}/COMPILED_OUTPUT/", mode: "copy"
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
