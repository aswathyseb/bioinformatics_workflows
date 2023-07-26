"""
This pipeline calls variants from the data.
"""

import csv, os, sys
from ruffus import *
import ruffus.cmdline as cmdline
import pandas as pd
import yaml, os
import helper

# Read config and get all parameters
args = helper.get_bwa_args()

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


if __name__ == "__main__":
    helper.run(args)
