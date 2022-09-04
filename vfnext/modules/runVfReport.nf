process runVfReport {
    publishDir "${params.outDir}/"

    input:
        path(fastp_htmls)
        path(snpeff_htmls)

    output:
        file("VF_REPORT/")
    script:
    """
    vfreports --glob-fastp  "*fastp.html" \
              --glob-snpeff "*snpEff*.html" \
              ./ ./VF_REPORT/
    """
}