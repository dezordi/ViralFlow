
process prepareDatabase {
    /*
    * Indexes reference fasta file using bwa.
    */
    publishDir "${params.databaseDir}/", mode:"copy"
    label "singlethread"
    input:
        val(genome_code)

    output:
        path("${genome_code}.fa"), emit: ref_fa
        path("${genome_code}.gff"), emit: ref_gff

    script:
        db_path = params.databaseDir
        base_path = "${db_path}/${genome_code}"
        // if files already available, just create symlinks
        if (file("${base_path}.fa").exists()){
          """
          echo \${PWD}/${genome_code}.fa
          ln -s ${base_path}.fa \${PWD}/${genome_code}.fa

          echo \${PWD}/${genome_code}.gff
          ln -s ${base_path}.gff \${PWD}/${genome_code}.gff
          """
        }
        // if gff or fa not available, download it
        else{
        """
          # download fasta from NCBI
          esearch -db nucleotide -query ${genome_code} \
             | efetch -format fasta > ${genome_code}.fa
          # download gff from NCBI
          wget -O ${genome_code}.gff \
           "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=${genome_code}"
        """
      }
}
