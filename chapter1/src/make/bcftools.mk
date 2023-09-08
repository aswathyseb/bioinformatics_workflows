
# The reference genome.
REF ?= refs/genome.fa

# The alignment file.
BAM ?= bam/S1.bam

# The variant file.
VCF ?= results/vcf/$(notdir $(basename ${BAM})).vcf.gz

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
usage:
> @echo "#"
> @echo "# Calls variants using bcftools"
> @echo "#"
>@echo "make vcf REF=${REF} BAM=${BAM}"
> @echo "#"
> @echo "# REF: Reference Genome"
> @echo "# BAM: Alignment file"
> @echo "#"

${VCF}: ${BAM} ${REF}
> mkdir -p $(dir $@)
>bcftools mpileup ${PFLAG} -O u -f ${REF} ${BAM} | \
 bcftools call ${CFLAG} -mv -O u | \
 bcftools norm -f ${REF} -d all -O u | bcftools sort -O z > ${VCF}

${VCF}.tbi: ${VCF}
> bcftools index -t -f ${VCF}

vcf: ${VCF}.tbi
> @ls -lh ${VCF}

all: vcf
> @ls -lh ${VCF}
