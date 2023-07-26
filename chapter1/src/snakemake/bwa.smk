import pandas as pd
import os

configfile: "config/config.yaml"


IDX_DIR = config.get('idx_dir')
ALN_DIR = config.get('aln_dir')
NCPU = config.get('ncpu')

BFLAGS = config.get('aligner').get('aln_flags')
SFLAGS = config.get('sam').get('sam_flags')
REF = config.get('genome')
IDX = f"{IDX_DIR}/{os.path.split(REF)[1]}"
IDX_FILE = f"{IDX}.bwt"

SHEET = config.get('data').get("sample_sheet")
LIB = config.get('data').get("library")

SAMPLES = pd.read_csv(SHEET,sep=",").set_index('sample',drop=False)


def get_reads_PE(wildcards):
    r1 = SAMPLES.loc[wildcards.sample, "read1"]
    r2 = SAMPLES.loc[wildcards.sample, "read2"]
    return r1, r2


def get_reads_SE(wildcards):
    return SAMPLES.loc[wildcards.sample, "read1"]


rule all_bams:
    input:
        IDX_FILE,
        expand(ALN_DIR + "/{sample}.bam",sample=SAMPLES.index),
        expand(ALN_DIR + "/{sample}.bam.bai",sample=SAMPLES.index)


rule index_genome:
    input:
        REF
    output:
        IDX_FILE
    shell:
        "bwa index -p {IDX} {input}"

if LIB == "SE":
    rule align_SE:
        input:
            IDX_FILE,
            reads=get_reads_SE
        output:
            ALN_DIR + "/{sample}.bam"
        threads: NCPU
        params:
            rg=r"@RG\tID:{sample}\tSM:{sample}"
        shell:
            "bwa mem -t {threads} {BFLAGS} -R '{params.rg}' {IDX} {input.reads} | samtools view -h {SFLAGS} |"
            "samtools sort >{output}"

if LIB == "PE":
    rule align_PE:
        input:
            IDX_FILE,
            reads=get_reads_PE
        output:
            ALN_DIR + "/{sample}.bam"
        threads: NCPU
        params:
            rg=r"@RG\tID:{sample}\tSM:{sample}"
        shell:
            "bwa mem -t {threads} {BFLAGS} -R '{params.rg}' {IDX} {input.reads} | samtools view -h {SFLAGS} |"
            "samtools sort >{output}"

rule index_bam:
    input:
        ALN_DIR + "/{sample}.bam"
    output:
        ALN_DIR + "/{sample}.bam.bai"
    shell:
        "samtools index {input}"
