# VFNext

Implementation of ViralFlow on Nextflow.

## current dev stage: alpha
A first version to run Viralflow under Nextflow framework was developed, tested on small dataset and everything works as expected.
Current development stage (alpha) requires investigation on real world scenarios to guarantee that VFNext implementation is equivalent to the first version of the pipeline.

## quick and dirty guide (only for dev purpose)

* How to setup vfnext?

```{bash}
git clone https://github.com/dezordi/ViralFlow.git
cd ViralFlow
git checkout vfnext
cd ./vfnext/containers/
/bin/bash setupContainers.sh
```

* How to run (on SARS-CoV-2)?

```{bash}
mkdir myRun
cd MyRun
nextflow run /path/to/vfnext/main.nf \
        --inDir /path/to/input_dir/ \
        --outDir /path/to/output_dir/ \
        --virus sars-cov2
        --adaptersFile /path/to/adapters.fa
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
        --adaptersFile /path/to/adapters.fa \
        --referenceGFF /path/to/reference.gff3
        --referenceGenome /path/to/reference.fasta
        --refGenomeCode my_genome_code
```

---
## NOTES
Paths provided for the parameters **must be absolute paths**

Pangolin, Nextclade and compile output steps will only run when virus tag is set to 'sars-cov2'

---
## Future Features
- [ ] Add an auto install procedure
- [ ] add CLI integration (similar to previous ViralFlow version)
- [ ] Add support for other virus (currently only 'sars-cov2' and 'custom' mode are supported)

## Roadmap
1. [x] ~Prototype~ (06/06/22)
2. [ ] Alpha
3. [ ] Beta
4. [ ] Release - ViralFlow v1.0
