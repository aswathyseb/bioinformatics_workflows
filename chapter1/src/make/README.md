### Description

This folder contains a variant calling workflow implemented as a Makefile.


### Workflow details

The main workflow is in `variants.mk` file. The workflow aligns data to the reference genome and calls variants.

The modules included in the variant caller workflow are `bwa.mk` and `bcftools.mk`.

`bwa.mk` : Alignment workflow for both single-end and paired-end reads using bwa-mem program.

`bcftools.mk` : Variant calling workflow using bcftools.

The master variant calling workflow `variants.mk` combines `bwa.mk` and `bcftools.mk`.

### How to run?

      make -f src/make.variants.mk REF=refs/genome.fa SAMPLE=S1 R1=data/S1_R1.fq R2=data/S1_R2.fq

The command will align sample `S1` in the data directory and call variants. By default the outputs will be in `results` folder.


User can specify the output bam and vcf file as below

   
      make -f src/make.variants.mk REF=refs/genome.fa SAMPLE=S1 R1=data/S1_R1.fq R2=data/S1_R2.fq BAM=S1.bam VCF=S1.vcf.gz


To call variants from multiple samples, combine the Makefile workflow with parallel.

      cat config/samples.csv | parallel --verbose --colsep , --header : make -f src/make/variants.mk  SAMPLE={sample} R1={read1} LIB=SE all


To combine the multiple vcf files into a single one, use the command

      make -f src/make/variants.mk VCF_TYPE=multi-sample MULTI_VCF=all.vcf.gz VCF_DIR=vcfs

To view the usage 
    
      make -f src/make/variants.mk

To perform a dry-run

    make -f src/make/variants.mk all --dry-run

### Sample sheet specification

Sample sheet is a comma-separated file with the header sample,read1,read2

An example sample sheet would look like this:

    sample,read1,read2
    S1,data/S1_R1.fq,data/S1_R2.fq
    S2,data/S2_R1.fq,data/S2_R2.fq

### Workflow requirements

1. Software requirements
   1. bwa
   2. bcftools
   3. gnu-parallel
2. Configuration requirements
   1. A samplesheet in comma-separated format with atleast 2 columns with column headers as 'sample,read1,read2'



