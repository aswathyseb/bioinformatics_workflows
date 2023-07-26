"""
This pipeline calls variants from the data.
"""

import csv, os, sys
from ruffus import *
import ruffus.cmdline as cmdline
import pandas as pd
import helper

# Read config and get all parameters
args = helper.get_variant_args()

# Get sample information
SAMPLES = pd.read_csv(args.samples, sep=",").set_index('sample', drop=False)

# Read data into a list of lists
reads = helper.read_data(SAMPLES, args.lib)

# Build the pattern from data directory and file suffix.
# All files should be in the same directory and should have the same suffix.
# Suffix is everything after the last underscore.
#
file1 = SAMPLES.iloc[0, 1]
sample1 = SAMPLES.index[0]
data_dir, s1 = helper.get_prefix_suffix(file1, sample1)
patt = f"{data_dir}/(.*){s1}"


# Task1 : Index reference genome
@transform(args.genome, regex(args.genome), args.idx_file)
def index_genome(inp, out):
    os.makedirs(os.path.dirname(out), exist_ok=True)
    cmd = f"bwa index -p {args.index} {inp}"
    print(cmd)
    os.system(cmd)


# Task2 : Align data to reference genome
@follows(index_genome)
@transform(reads, regex(patt), add_inputs(args.genome), args.aln_dir + "/\\1.bam")
def align_reads(inp, out):
    os.makedirs(args.aln_dir, exist_ok=True)
    file1 = inp[0][0]
    if args.lib == "PE":
        # Align in paired end mode
        file2 = inp[0][1]
        cmd = f"bwa mem -t {args.jobs} {args.aln_flags} {args.index} {file1} {file2} |samtools view -h {args.sam_flags} |samtools sort > {out}"
    if args.lib == "SE":
        # Align in single end mode
        cmd = f"bwa mem -t {args.jobs}  {args.aln_flags} {args.index} {file1}  | samtools view -h {args.sam_flags} | samtools sort > {out}"
    print(cmd)
    os.system(cmd)


# Task3 : Index bam files.
@follows(align_reads)
@transform(align_reads, suffix(".bam"), ".bam.bai")
def index_bam(inp, out):
    cmd = f"samtools index {inp}"
    print(cmd)
    os.system(cmd)


MULTI_VCF = False
SAMPLE_VCF = False
if args.vcf_type == "multi-sample":
    MULTI_VCF = True
if args.vcf_type == "sample-vcf":
    SAMPLE_VCF = True


# Task4 : Call variants for each sample
@active_if(SAMPLE_VCF)
@follows(align_reads)
@transform(align_reads, formatter(), args.vcf_dir + "/{basename[0]}.vcf.gz")
def create_sample_vcf(inp, out):
    os.makedirs(args.vcf_dir, exist_ok=True)
    cmd = f"""
            bcftools mpileup {args.pflags} -O u -f {args.genome} {inp} |
            bcftools call {args.cflags} -mv -O u |
            bcftools norm -f {args.genome} -d all -O u | bcftools sort -O z > {out}
            """
    print(cmd)
    os.system(cmd)


# Task 5 : Index vcf file
@active_if(SAMPLE_VCF)
@follows(create_sample_vcf)
@transform(create_sample_vcf, suffix(".vcf.gz"), ".vcf.gz.csi")
def index_sample_vcf(inp, out):
    cmd = f"""bcftools index {inp}"""
    print(cmd)
    os.system(cmd)


# Task4 : Create a multi-sample vcf
@active_if(MULTI_VCF)
@follows(align_reads)
@merge(align_reads, args.vcf_dir + "/all.vcf.gz")
def create_multisample_vcf(inp, out):
    os.makedirs(args.vcf_dir, exist_ok=True)
    inputs = " ".join(inp)
    cmd = f"""bcftools mpileup {args.pflags} -O u -f {args.genome} {inputs} | 
               bcftools call {args.cflags} -mv -O u | 
               bcftools norm -f {args.genome} -d all -O u | bcftools sort -O z > {out}
            """
    print(cmd)
    os.system(cmd)


# Task 5 : Index vcf file
@active_if(MULTI_VCF)
@follows(create_multisample_vcf)
@transform(create_multisample_vcf, suffix(".vcf.gz"), ".vcf.gz.csi")
def index_multisample_vcf(inp, out):
    cmd = f"""bcftools index {inp}"""
    print(cmd)
    os.system(cmd)


if __name__ == "__main__":
    helper.run(args)
