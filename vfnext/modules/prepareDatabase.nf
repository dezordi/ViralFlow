
process prepareDatabase {
    /*
    * Indexes reference fasta file using bwa.
    */
    publishDir "${params.databaseDir}/"

    input:
        val(genome_code)

    output:
        path("*.fa"), emit: ref_fa
        path("*.gff"), emit: ref_gff

    script:
        db_path = params.databaseDir
        """
        # check if files are available at database dir,
        # if not download it
        if [ -f "${db_path}/${genome_code}.fa" ];then
          echo "${db_path}/${genome_code}.fa alread exists."
          ln -s ${db_path}/${genome_code}.fa ./
        else
          esearch -db nucleotide -query ${genome_code} \
             | efetch -format fasta > ./${genome_code}.fa
        fi

        if [ -f "${db_path}/${genome_code}.gff" ]; then
          echo "${db_path}/${genome_code}.gff alread exists."
          ln -s ${db_path}/${genome_code}.gff ./

        else
          wget -O ./${genome_code}.gff \
           "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=${genome_code}"
        fi
        """
}
