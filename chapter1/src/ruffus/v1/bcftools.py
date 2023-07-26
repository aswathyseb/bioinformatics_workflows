"""
This pipeline calls variants from the data.
"""

import csv, os, sys
from ruffus import *
import ruffus.cmdline as cmdline
import pandas as pd
import yaml, os, glob
import helper


def bcftools_pipeline(args):
    SAMPLES = []
    if args.samples:
        # Get sample information
        SAMPLES = pd.read_csv(args.samples, sep=",").set_index('sample', drop=False)
        SAMPLES = SAMPLES.index.to_list()

    bams = helper.read_bam(args.aln_dir, samples=SAMPLES)

    if not bams:
        print("No input files given!!")
        sys.exit()

    MULTI_VCF = False
    SAMPLE_VCF = False
    if args.vcf_type == "multi-sample":
        MULTI_VCF = True
    if args.vcf_type == "sample-vcf":
        SAMPLE_VCF = True

    # Task4 : Call variants for each sample
    @active_if(SAMPLE_VCF)
    @transform(bams, formatter(), args.vcf_dir + "/{basename[0]}.vcf.gz")
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
    @merge(bams, args.vcf_dir + "/all.vcf.gz")
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
    # Read config and get all parameters
    args = helper.get_bcftools_args()
    bcftools_pipeline(args)
    helper.run(args)
