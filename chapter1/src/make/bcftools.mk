
# The reference genome.
REF ?= refs/genome.fa

# Alignment folder.
#ALIGN_DIR ?= results/bam

# Sample name
SAMPLE ?= S1

# The alignment file.
BAM ?= results/bam/${SAMPLE}.bam

# VCf file
VCF ?= results/vcf/${SAMPLE}.vcf.gz

# Pileup flags
PFLAG ?= -d 100

# Variant calling flags
CFLAG ?= --ploidy 2

# Makefile customizations.
.RECIPEPREFIX = >
.DELETE_ON_ERROR:
.ONESHELL:
MAKEFLAGS += --warn-undefined-variables --no-print-directory

#Print usage.
bcftools_usage:
> @echo "#"
> @echo "# Calls variants using bcftools"
> @echo "#"
>@echo "make vcf REF=${REF} SAMPLE=${SAMPLE} BAM=${BAM} VCF=vcf/${SAMPLE}.vcf.gz"
> @echo "#"
> @echo "# REF       : Reference Genome"
> @echo "# SAMPLE    : Sample name"
> @echo "# ALIGN_DIR : Alignment file"
> @echo "# VCF       : Path to output vcf"
> @echo "#"


${VCF}: ${BAM}.bai ${REF}
> echo ${BAM}
> mkdir -p $(dir $@)
>bcftools mpileup ${PFLAG} -O u -f ${REF} ${BAM} | \
 bcftools call ${CFLAG} -mv -O u | \
 bcftools norm -f ${REF} -d all -O u | bcftools sort -O z > ${VCF}

${VCF}.tbi: ${VCF}
> bcftools index -t -f ${VCF}

vcf: ${VCF}.tbi
> @ls -lh ${VCF}
