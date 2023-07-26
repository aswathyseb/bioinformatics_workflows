### Description

This folder contains snakemake variant calling workflow.


### Workflow details

The main workflow is in `variants.smk` file. The workflow aligns data to the reference genome and calls variants.

The modules included in the variant caller workflow are `bwa.smk` and `bcftools.smk`.

`bwa.smk` : Alignment workflow for both single-end and paired-end reads using bwa-mem program.

`bcftools.smk` : Variant calling workflow using bcftools.

The master variant calling workflow `variants.smk` combines `bwa.smk` and `bcftools.smk`.

### How to run?

Run using 4 processors and parameters specified in config file.

      snakemake -s src/snakemake/variants.smk -c 4

To override any of the runtime parameters specified in the config file, for eg: alignment directory

      snakemake -s src/snakemake/variants.smk -c 4 --config aln_dir=align

To create a multi-sample vcf, run the command

      snakemake -s src/snakemake/variants.smk -c 4 --config vcf_type=multi-sample
   
To call variants from each sample and create a sample specific vcf, run the command

      snakemake -s src/snakemake/variants.smk -c 4 --config vcf_type=sample-vcf

To run with a different config file

      snakemake -s src/snakemake/variants.smk -c 4 --configfile new_config.yaml

To view all the commandline options
    
      snakemake -h

To perform a dry-run

    python src/ruffus/vcf.py -n 

To perform a dry-run

    snakemake -s src/snakemake/variants.smk -c 4 --config aln_dir=align -n

### Sample sheet specification

Sample sheet is a comma-separated file with the header sample,read1,read2

An example sample sheet would look like this:

    sample,read1,read2
    S1,data/S1_R1.fq,data/S1_R2.fq
    S2,data/S2_R1.fq,data/S2_R2.fq

### Workflow requirements

1. Software requirements 
   1. Snakemake
   2. bwa
   3. bcftools
2. Configuration requirements
   1. All parameters are specified through a YAML configuration file, eg: config.yaml 
   2. A samplesheet in comma-separated format with atleast 2 columns with column headers as 'sample,read1,read2'



