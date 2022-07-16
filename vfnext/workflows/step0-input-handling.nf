#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {prepareDatabase} from "../modules/prepareDatabase.nf"

// set supported virus flag

def validate_parameters() {
    // --- SANITY CHECKS ------------------------------------------------------
    def errors = 0
    // check if required params were provided and if files provided exists

    if (params.adaptersFile==null){
      // TODO: make adapter file optional, usefull for metagenomics
      log.error("An adapter file path must be provided")
      errors +=1
    }
    // --- VIRUS FLAGS CHECK
    // check if a valid virus flag was provided
    def valid_virus = ["sars-cov2","custom"]
    if (!valid_virus.contains(params.virus)) {
      log.error("The virus provided (${params.virus}) is not valid.")
      errors += 1
    }

    // be sure custom only options were not set if a valid virus tag was provided
    if (valid_virus.contains(params.virus) && !(params.virus == "custom")) {
        if (!(params.referenceGFF==null)){
          log.warn("The valid virus tag (${params.virus}) was provided, ingnoring the provided referenceGFF (${params.referenceGFF})")
          params.referenceGFF=null
        }
        if (!(params.referenceGenome==null)){
          log.warn("The valid virus tag (${params.virus}) was provided, ingnoring the provided referenceGenome (${params.referenceGenome})")
          params.referenceGenome=null
        }
        if (!(params.refGenomeCode==null)){
          log.warn("The valid virus tag (${params.virus}) was provided, ingnoring the provided refGenomeCode (${params.refGenomeCode})")
          params.refGenomeCode=null
        }

    }
    // ------------------------------------------------------------------------
    // if a custom virus, check if mandatory params were set
    if (params.virus=="custom"){
      // if a genome code was not provided, check if a gff and a ref fasta was
      if (params.refGenomeCode==null){
        if (params.runSnpEff==true){
          log.warn("The runSnpEff was set to ${params.runSnpEff}, but no refGenomeCode was provided.")
          log.warn("SnpEff will not be run")
        }
        if (params.referenceGFF==null){
          log.error("A 'custom' virus tag was set and no refGenomeCode was provided, therefore a referenceGFF must be provided.")
          errors += 1
        } else {
          ref_gff_path = file(params.referenceGFF)
          if (!ref_gff_path.isFile()){
            log.error("${ref_gff_path} is not a file.")
            errors += 1
          }
          if (!ref_gff_path.exists()){
            log.error("${ref_gff_path} does not exists.")
            errors += 1
          }
        }

        if (params.referenceGenome==null){
          log.error("A 'custom' virus tag was set and no refGenomeCode was provided, therefore a referenceGenome must be provided.")
          errors += 1
        } else {
          ref_fa_path = file(params.referenceGenome)
          if (!ref_fa_path.isFile()){
            log.error("${ref_path} is not a file.")
            errors += 1
          }
          if (!ref_fa_path.exists()){
            log.error("${ref_path} does not exists.")
            errors += 1
          }
        }
      }
    }
    // ------------------------------------------------------------------------

    // check if output dir exists, if not create the default
    if (params.outDir){
       outDir_path = file(params.outDir)

       if (!outDir_path.isDirectory()){
         log.error("${params.outDir} is not a directory")
         error+=1
       }

       if (!outDir_path.exists()){
         log.warn("${params.outDir} does not exist, the directory will be created")
         outDir_path.mkdir()
       }
    }
    // check if input dir exists
    if (params.inDir==null){
        log.error("An input directory must be provided.")
        error+=1
    }
    if (params.inDir){
      inDir_path = file(params.inDir)
      if (!inDir_path.isDirectory()){
        log.error("${params.inDir} is not a directory")
        error+=1
      }

    }

    //TODO check if adapterFile exist
    //-------------------------------------------------------------------------
    // count errors and kill nextflow if any had been found
    if (errors > 0) {
        log.error(String.format("%d errors detected", errors))
        exit 1
    }
}

workflow processInputs {
  main:
    // --- Sanity Check -------------------------------------------------------
    // check if fasta exists and follow symlinks if needed
    //ref_fa = file(reference_fasta, checkIfExists=true, followLinks=true)
    //-------------------------------------------------------------------------
    validate_parameters()

    // ---- get reference GFF and fasta ---------------------------------------
    // Setup ref code values for supported virus
    ref_gcode = null
    reference_fa = null
    reference_gff = null

    if (!(params.virus=="custom")){
      if (params.virus=="sars-cov2"){
        ref_gcode = "NC_045512.2"
      }
    }

    // if custom virus, check if a genome code was provided, if not
    // emit the ref gff and fasta provided
    if (params.virus=="custom"){
      if (!(params.refGenomeCode==null)){
        ref_gcode = params.refGenomeCode
      } else {
        reference_gff = params.referenceGFF
        reference_fa = params.referenceGenome
      }
    }

    // if a genome code was provided, get the reference fasta and gff
    if (!(ref_gcode==null)){
      prepareDatabase(ref_gcode)
      reference_fa = prepareDatabase.out.ref_fa
      reference_gff = prepareDatabase.out.ref_gff
    }

    // get reads
    reads_channel_fqgz = channel.fromFilePairs("${params.inDir}/*_R{1,2}*.fq.gz")
    reads_channel_fagz = channel.fromFilePairs("${params.inDir}/*_R{1,2}*.fastq.gz")
    reads_channel = reads_channel_fagz.concat(reads_channel_fqgz)

    if (reads_channel.count()==0){
      log.error("No fastq files found. Be sure your fastq can be found by '*_R{1,2}*.fq.gz'(or fasta.gz) wildcard.")
      exit 1
    }
    // be sure a reference fasta and a reference gff was obtained
    assert !(reference_fa == null) && !(reference_gff == null)

  emit:
    reads_ch = reads_channel
    ref_gff = reference_gff
    ref_fa = reference_fa
    ref_gcode = ref_gcode
}
