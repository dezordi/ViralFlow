#!/usr/bin/env nextflow


// -----------------------------------------------------------------------------
// TODO: - create container recipe for ivar + samtools
//       - Check with Filipe the if ivar step can be done just after the bwa index
// -----------------------------------------------------------------------------

// enable dsl2
nextflow.enable.dsl = 2

// import modules
include { indexReferenceBWA } from './modules/bwaIndex.nf'
include { runFastp } from './modules/fastp.nf'
include { align2ref } from './modules/align2ref.nf'
include { runIvar } from './modules/runIvar.nf'

// I got some of the code from the FASTQC PIPELINE (https://github.com/angelovangel/nxf-fastqc/blob/master/main.nf)

/*
* ANSI escape codes to color output messages, get date to use in results folder name
*/
ANSI_GREEN = "\033[1;32m"
ANSI_RED = "\033[1;31m"
ANSI_RESET = "\033[0m"

log.info """
        ===========================================
         VFNEXT_0.0 (dev : pre-alpha)
         Used parameters:
        -------------------------------------------
         --inDir            : ${params.inDir}
         --outDir           : ${params.outDir}
         --referenceGenome  : ${params.referenceGenome}
         --adaptersFile     : ${params.adaptersFile}
         --threads          : ${params.threads}
         --minLen           : ${params.minLen}
         --depth            : ${params.depth}
         --minDpIntrahost   : ${params.minDpIntrahost}
         --trimLen          : ${params.trimLen}

         Runtime data:
        -------------------------------------------
         Running with profile:   ${ANSI_GREEN}${workflow.profile}${ANSI_RESET}
         Run container:          ${ANSI_GREEN}${workflow.container}${ANSI_RESET}
         Running as user:        ${ANSI_GREEN}${workflow.userName}${ANSI_RESET}
         Launch dir:             ${ANSI_GREEN}${workflow.launchDir}${ANSI_RESET}
         Base dir:               ${ANSI_GREEN}${baseDir}${ANSI_RESET}
         ------------------------------------------
         """
         .stripIndent()


//  The default workflow
workflow {
   //println "\nI want to do the genome indexing of $params.referenceGenome and put the output at $params.outDir"
   // STEP 1 ------------------------------------------------------------------
   // open input channels
   // for bwa index
   refgen_ch = channel.fromPath(params.referenceGenome)
   //refgen_ch.view()
   reads_ch = channel.fromFilePairs('input/*_R{1,2}*.fq.gz')
   //reads_ch.view()

   // run indexing, open the bwa index output channel
   indexReferenceBWA(refgen_ch)
   indexReferenceBWA.out.set { bwaidx_Output_ch }
   //bwaidx_Output_ch.view()
   // run fastp
   runFastp(reads_ch)
   runFastp.out.set { fastp_Output_ch }
   //fastp_Output_ch.view()

   //align 2 reference
   align2ref_In_ch = reads_ch.combine(bwaidx_Output_ch)
   align2ref(align2ref_In_ch)
   align2ref.out.set { align2ref_Out_ch }
   // ivar
   ivar_In_ch = channel.from(params.referenceGenome, params.depth, align2ref_Out_ch)
   runIvar(align2ref_Out_ch)
}

// -------------- Check if everything went okay -------------------------------
workflow.onComplete {
    if (workflow.success) {
        log.info """
            ===========================================
            ${ANSI_GREEN}Finished in ${workflow.duration}
            See the report here ==> ${ANSI_RESET}$params.outDir/XXX_report.html
            """
            .stripIndent()
    } else {
        log.info """
            ===========================================
            ${ANSI_RED}Finished with errors!${ANSI_RESET}
            """
            .stripIndent()
    }
}
