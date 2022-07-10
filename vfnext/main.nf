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
include { runReadCounts } from './modules/runReadCounts.nf'
include { alignConsensus2Ref } from './modules/alignConsensus2Ref.nf'
include { runIntraHostScript } from './modules/runIntraHostScript.nf'
include { runPangolin } from './modules/runPangolin.nf'
include { runNextClade } from './modules/runNextclade.nf'
include { runPicard } from './modules/runPicard.nf'
include { fixWGS } from './modules/fixWGS.nf'
include { compileOutputs } from './modules/compileOutput.nf'
include { processInputs } from './workflows/step0-input-handling.nf'
include { genFaIdx } from './modules/genFaIdx.nf'

// I got some of the code from the FASTQC PIPELINE
// https://github.com/angelovangel/nxf-fastqc/blob/master/main.nf

/*
* ANSI escape codes to color output messages, get date to use in results folder name
*/
ANSI_GREEN = "\033[1;32m"
ANSI_RED = "\033[1;31m"
ANSI_RESET = "\033[0m"

log.info """
        ===========================================
         VFNEXT_0.1 (dev : alpha)
         Used parameters:
        -------------------------------------------
         --inDir            : ${params.inDir}
         --outDir           : ${params.outDir}
         --virus            : ${params.virus}
         --referenceGenome  : ${params.referenceGenome}
         --referenceGFF     : ${params.referenceGFF}
         --adaptersFile     : ${params.adaptersFile}
         --threads          : ${params.threads}
         --minLen           : ${params.minLen}
         --depth            : ${params.depth}
         --minDpIntrahost   : ${params.minDpIntrahost}
         --trimLen          : ${params.trimLen}
         --databaseDir      : ${params.databaseDir}
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
   // STEP 0 ------------------------------------------------------------------
   // open input channels
   processInputs()
   reads_ch = processInputs.out.reads_ch
   ref_gff = processInputs.out.ref_gff
   ref_fa = processInputs.out.ref_fa
   ref_gcode = processInputs.out.ref_gcode

   // STEP 1 ------------------------------------------------------------------
   // run indexing, open the bwa index output channel
   indexReferenceBWA(ref_fa)
   indexReferenceBWA.out.set { bwaidx_Output_ch }
   genFaIdx(ref_fa)
   genFaIdx.out.set {faIdx_ch}

   // run fastp
   runFastp(reads_ch)

   //align 2 reference
   align2ref_In_ch = reads_ch.combine(bwaidx_Output_ch)
   align2ref(align2ref_In_ch, ref_fa)
   align2ref.out.set { align2ref_Out_ch }
   // ivar
   runIvar(align2ref_Out_ch, ref_fa)
   runIvar.out.set { runIvar_Out_ch }

   // readcounts
   runReadCounts(align2ref_Out_ch)
   runReadCounts.out.set {runReadCounts_Out_ch}

   //align consensus to ref
   alignConsensus2Ref(runIvar_Out_ch, ref_fa)
   alignConsensus2Ref.out.set {alignCon_Out_ch}

   // Assembly Metrics
   runPicard(align2ref_Out_ch, ref_fa)
   runPicard.out.set {runPicard_Out_ch}
   fixWGS_In_ch =runPicard_Out_ch.join(runIvar_Out_ch)
   fixWGS(fixWGS_In_ch)

   //run intrahost
   intraHost_In_ch = alignCon_Out_ch.join(runReadCounts_Out_ch)
   runIntraHostScript(intraHost_In_ch, ref_gff)
   runIntraHostScript.out.set {runIntraHostScript_Out_ch}

   // run Variant Naming (Pangolin and Nextclade)
   runVariantNaming_In_ch = runIntraHostScript_Out_ch.join(runIvar_Out_ch)
   runPangolin(runVariantNaming_In_ch)
   runNextClade(runVariantNaming_In_ch, ref_fa)

   // GAMBIARRA ALLERT --------------------------------------------------------
   // Pangolin is the last ones to run, so will use it as a trigger to
   // the output compilation/
   // for the final version, need to find a better way. Maybe split and set
   // as individual post analysis workflow
   final_trigger = runPangolin.out.collect()
   compileOutputs(final_trigger)

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
