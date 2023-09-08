
#
# Call variants from data based on a reference
#

include src/make/bwa.mk
include src/make/bcftools.mk

# Type of variant calling multi-sample or per sample vcf.
# Options are multi-sample or sample-vcf
VCF_TYPE ?= multi-sample

VCF_FILE ?= results/vcf/all.vcf.gz

ALIGN_DIR = results/bam

# Makefile customizations.
.RECIPEPREFIX = >
.DELETE_ON_ERROR:
.ONESHELL:
MAKEFLAGS += --warn-undefined-variables --no-print-directory


#Print usage.
usage:
> @echo "#"
> @echo "# Calls variants from data based on a reference"
> @echo "#"
>@echo "make vcf REF=${REF} SAMPLE=${SAMPLE} LIB=${LIB} R1=${R1} R2=${R2} VCF_TYPE=${VCF_TYPE}"
> @echo "#"
> @echo "# REF    : Reference Genome"
> @echo "# LIB    : Library type, PE or SE"
> @echo "# SAMPLE : Sample name"
> @echo "# R1     : Read1"
> @echo "# R2     : Read2"
> @echo "#"

