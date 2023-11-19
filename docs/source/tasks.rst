Tasks
=====

.. _prepare_inputs:

Prepare Inputs
--------------

In this step the paired reads are identified on input directory, and if the user provides a NCBI genome code, the fasta and gff file of the respective genome will be downloaded.


.. _quality_control:

Quality Control
---------------

In this step, the paired reads are processed using the `fastp <https://github.com/OpenGene/fastp>`_, with the following parameters:

* **--detect_adapter_for_pe**: To enable auto detection of adapters.
* **-l <minLen>**: To remove reads with length lower than the threshold parsed on `minLen` argument.
* **-f -t -F -T <trimLen>**: To remove reads boundaries regions, based on `trimLen`.
* **--cut_front --cut_tail --qualified_quality_phred 20**: Drop bases in tail/fron based on a sliding windown (4 bases) with mean quality less than 20.

.. _mapping_to_reference:

Mapping to Reference
--------------------

In this step, the processed reads in `Quality Control` step are mapped against the reference genome using `BWA <https://bio-bwa.sourceforge.net/>`_, with default parameters. If a bed files with PCR primers positions is parsed, the primers regions are clipped using `Samtools <https://github.com/samtools/samtools>`_ ampliconclip module. 

.. _get_map_unmapped_reads:

Get Map/Unmapping Reads
-----------------------

(Optional) In this step, the mapped and/or unmapped reads are extracted to specific fastq files using `Samtools <https://github.com/samtools/samtools>`_.

.. _generate_consensus:

Generate Consensus
------------------

In this step, the consensus files are generated with `Samtools <https://github.com/samtools/samtools>`_ and `iVar <https://github.com/andersen-lab/ivar>`_, with the following parameters:

* **-q <mapping_quality>**: Minimum quality score threshold to count base.
* **-t 0/0.6**: Minimum frequency threshold, one consensus is generated with major allele frequencies (-t 0) and a ambiguoues consensus is generated with IUPAC Ambiguity Codes for loci with 2 alleles in ~0.6 and ~0.4 alleles.
* **-m <depth>**: Minimum depth to call consensus, loci without the minimum depth is hard-masked with N.

.. _coverage_plot:

Coverage Plot
-------------

In this step, a genomic map showing the coverage depth is generated using the sorted BAM file produced in the `Mapping to Reference` step. The tool `BAMdash` is used for this purpose.

.. _snp_plot:

SNP Plot
--------

In this step, a plot is generated to visualize the positions of the variants found in samples. This provides an overview of the distribution of single nucleotide polymorphisms (SNPs) across the genome. The tool `snipit` is used for this purpose.

.. _variant_annotation:

Variant Annotation
------------------

(Optional) In this step, the variants are called with freebayes `Freebayes <https://github.com/freebayes/freebayes>`_ (default parameters) and annotated with `snpEff <https://pcingola.github.io/SnpEff/>`_. 

.. _get_metrics:

Get Assembly Metrics
--------------------

In this step, the assembly metrics are evaluated using `picard <https://broadinstitute.github.io/picard/>`_.

.. _get_intrahosts:

Get Intrahosts
--------------

In this step, the `bam-readcount <https://github.com/genome/bam-readcount>`_. and an `in house` python script are used to identify intrahosts regions based on different rules:

* A minimum (minDpIntrahost) depth do consider intrahosts regions
* The intrahost in minor allele frequency should represent at least 5% of loci total depth.
* The intrahost in minor allele frequency should have at least 5% in each read sense (foward and reverse)

.. _run_pangolin:

Run Pangolin
------------

(Optional) For SARS-CoV-2 analysis, the `Pangolin <https://github.com/cov-lineages/pangolin>`_  tool is used to perform the PANGO lineage signature.

.. _run_nextclade:

Run Nextclade
-------------

(Optional) For SARS-CoV-2 analysis, the `Nextstrain <https://github.com/nextstrain/nextclade>`_ tool is used to perform the Nexstrain clade signature, and generate different genome metrics.

.. _vf_report:

Get ViralFlow report
--------------------

In this step, an HTML report is created with fastp and snpEff results with  `ViralFlow-Report <https://github.com/dezordi/viralflow-report>`_ script.
