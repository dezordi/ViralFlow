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

    // be sure custom options were not set if a valid virus tag was provided
    if (valid_virus.contains(params.virus)) {
        if (!(params.referenceGFF==null)){
          log.warn("The valid virus tag (${params.virus}) was provided, ingnoring the provided referenceGFF (${params.referenceGFF})")
          params.referenceGFF=null
        }
        if (!(params.referenceGenome==null)){
          log.warn("The valid virus tag (${params.virus}) was provided, ingnoring the provided referenceGenome (${params.referenceGenome})")
          params.referenceGenome=null
        }
    }

    // if a custom virus, check if mandatory params were set
    if (params.virus=="custom"){

      if (params.refGenomeCode==null){
        log.error("For a 'custom' virus, a genomeCode must be provided.")
        errors += 1
      }

      if (params.referenceGFF==null){
        log.error("A 'custom' virus tag was set, a referenceGFF must be provided.")
        errors += 1
      }

      if (params.referenceGenome==null){
        log.error("For a 'custom' virus, a reference fasta must be provided.")
        errors += 1
      }

      else {
        ref_path = file(params.referenceGenome)
        if (!ref_path.isFile()){
          log.error("${ref_path} is not a file.")
          errors += 1
        }
        if (!ref_path.exists()){
          log.error("${ref_path} does not exists.")
          errors += 1
        }
      }
    }

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
    // Setup values for supported virus
    if (!(params.virus=="custom")){
      if (params.virus=="sars-cov2"){
      ref_genome_code = "NC_045512.2"
      }
      prepareDatabase(ref_genome_code)
      reference_fa = prepareDatabase.out.ref_fa
      reference_gff = prepareDatabase.out.ref_gff
    }

    // if custom virus, just emit the ref gff and fasta provided
    if (params.virus=="custom"){
      reference_gff = params.referenceGFF
      reference_fa = params.referenceGenome
    }

    // get reads
    reads_channel = channel.fromFilePairs("${params.inDir}/*_R{1,2}*.fq.gz")

  emit:
    reads_ch = reads_channel
    ref_gff = reference_gff
    ref_fa = reference_fa
}
