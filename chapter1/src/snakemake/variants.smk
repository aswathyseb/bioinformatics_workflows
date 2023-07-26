import pandas as pd

# Workflow parameters

configfile: "config/config.yaml"

VCF_DIR = config.get('vcf_dir')
VCF_CALL = config.get("vcf_type")
SHEET = config.get('data').get("sample_sheet")
SAMPLES = pd.read_csv(SHEET,sep=",").set_index('sample',drop=False)

# Load rules

include: "bwa.smk"
include: "bcftools.smk"


if VCF_CALL == "multi-sample":
    rule all:
        input:
            VCF_DIR + "/all.vcf.gz",
            VCF_DIR + "/all.vcf.gz.csi",

if VCF_CALL == "sample-vcf":
    rule all:
        input:
            expand(VCF_DIR + "/{sample}.vcf.gz",sample=SAMPLES.index),
            expand(VCF_DIR + "/{sample}.vcf.gz.csi",sample=SAMPLES.index)
