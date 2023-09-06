# VFNext

Implementation of ViralFlow on Nextflow.

## quick start guide

* How to setup vfnext?

```{bash}
git clone https://github.com/dezordi/ViralFlow.git
cd ViralFlow
cd ./vfnext/containers/
/bin/bash setupContainers.sh
/bin/bash add_entries_SnpeffDB.sh 
```

* How to run (on SARS-CoV-2)?

```{bash}
mkdir myRun
cd MyRun
nextflow run /path/to/vfnext/main.nf \
        --inDir /path/to/input_dir/ \
        --outDir /path/to/output_dir/ \
        --virus sars-cov2
        --primersBED /path/to/adapters.fa
```

* How to run on a arbitrary virus?
To run viralflow on a non-supported virus, user must provide:
1. a reference gff file
2. a reference fasta file
3. a genome code (if user wants to use snpEff)

```{bash}
nextflow run ../vfnext/main.nf \
        --inDir /path/to/input_dir/ \
        --outDir /path/to/output_dir/ \
        --virus custom \
        --primersBED /path/to/adapters.fa \
        --referenceGFF /path/to/reference.gff3
        --referenceGenome /path/to/reference.fasta
        --refGenomeCode my_genome_code
```

---
## NOTES
Paths provided for the parameters **must be absolute paths**

Pangolin, Nextclade and compile output steps will only run when virus tag is set to 'sars-cov2'

---

