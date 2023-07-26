import pandas as pd


configfile: "config/config.yaml"

ALN_DIR = config.get('aln_dir')
VCF_DIR = config.get('vcf_dir')
REF = config.get('genome')
VCF_CALL = config.get("vcf_type")

SHEET = config.get('data').get("sample_sheet")

SAMPLES = pd.read_csv(SHEET,sep=",").set_index('sample',drop=False)

if VCF_CALL == "multi-sample":
    rule call_variants:
        input:
            VCF_DIR + "/all.vcf.gz",
            VCF_DIR + "/all.vcf.gz.csi",

if VCF_CALL == "sample-vcf":
    rule call_variants:
        input:
            expand(VCF_DIR + "/{sample}.vcf.gz",sample=SAMPLES.index),
            expand(VCF_DIR + "/{sample}.vcf.gz.csi",sample=SAMPLES.index)


rule call_multi_sample_variants:
    input:
        fa=REF,
        bam=expand(ALN_DIR + "/{sample}.bam",sample=SAMPLES.index),
        bai=expand(ALN_DIR + "/{sample}.bam.bai",sample=SAMPLES.index)
    output:
        VCF_DIR + "/all.vcf.gz"
    params:
        pflags="-d 100",
        cflags="--ploidy 2"
    shell:
        "bcftools mpileup {params.pflags} -O u -f {input.fa} {input.bam} | "
        "bcftools call {params.cflags} -mv -O u | "
        "bcftools norm -f {input.fa} -d all -O u | bcftools sort -O z > {output}"


rule call_sample_variants:
    input:
        fa=REF,
        bam=ALN_DIR + "/{sample}.bam",
        bai=ALN_DIR + "/{sample}.bam.bai"
    output:
        VCF_DIR + "/{sample}.vcf.gz"
    params:
        pflags=config.get('vcf')['pileup_flags'],
        cflags=config.get('vcf')['call_flags']
    shell:
        "bcftools mpileup {params.pflags} -O u -f {input.fa} {input.bam} | "
        "bcftools call {params.cflags} -mv -O u | "
        "bcftools norm -f {input.fa} -d all -O u | bcftools sort -O z > {output}"


rule index_multi_sample_vcf:
    input:
        VCF_DIR + "/all.vcf.gz"
    output:
        VCF_DIR + "/all.vcf.gz.csi"
    shell:
        "bcftools index {input}"

rule index_sample_vcf:
    input:
        VCF_DIR + "/{sample}.vcf.gz"
    output:
        VCF_DIR + "/{sample}.vcf.gz.csi"
    shell:
        "bcftools index {input}"
