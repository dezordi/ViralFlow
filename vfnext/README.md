# VFNext

Implementation of ViralFlow on Nextflow framework.

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

* How to run?

```{bash}
mkdir myRun
cd MyRun
nextflow run /path/to/vfnext/main.nf \
        --inDir /path/to/input_dir/ \
        --outDir /path/to/output_dir/ \
        --referenceGenome /path/to/reference_genome.fasta \
        --adaptersFile /path/to/adapters.fa
```

## Roadmap
1. [x] ~Prototype~ (06/06/22)
2. [ ] Alpha
3. [ ] Beta
4. [ ] Release - ViralFlow v1.0
