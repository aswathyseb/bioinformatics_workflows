
#
# Call variants from data based on a reference
#


# Makefile customizations
.RECIPEPREFIX = >
.DELETE_ON_ERROR:
.ONESHELL:
MAKEFLAGS += --warn-undefined-variables --no-print-directory

#Print usage.
usage:
> @echo "#"
> @echo "# Calls variants from data based on a reference"
> @echo "#"
> @echo "make all REF=${REF} SAMPLE=${SAMPLE} LIB=${LIB} R1=${R1} R2=${R2} BAM=${BAM} VCF=${VCF}"
> @echo "#"
> @echo "# REF    	: Reference Genome"
> @echo "# LIB    	: Library type, PE or SE"
> @echo "# SAMPLE 	: Sample name"
> @echo "# R1     	: Read1"
> @echo "# R2     	: Read2"
> @echo "#"
> @echo "# Other variables are"
> @echo "#"
> @echo "# BAM		: Path to bam file."
> @echo "# VCF    	: Path to vcf file."
> @echo "# MULTI_VCF	: Path to multi-vcf file."
> @echo "# VCF_DIR	: Directory with multiple vcf files."
> @echo "# VCF_TYPE	: multi-sample if multi-sample vcf is needed"
> @echo "#"


# Include bwa.mk and bcftools.mk
include src/make/bwa.mk
include src/make/bcftools.mk

# Variant calling settings
VCF_TYPE ?= sample-vcf # or 'multi-sample'
MULTI_VCF ?= results/vcf/all.vcf.gz
VCF_DIR ?= $(dir ${VCF})

# All vcf files for merging.
ALL_VCFS := $(wildcard ${VCF_DIR}/*.vcf.gz)

# Number of threads to use
T ?= 4

# Rule for multi-sample VCF generation
$(MULTI_VCF): $(wildcard ${VCF_DIR}/*.vcf.gz.tbi)
> mkdir -p $(dir $@)
> bcftools merge --threads $T -O z -o $@ $(ALL_VCFS)
> bcftools index -t -f ${MULTI_VCF}

# Conditional execution based on VCF_TYPE
ifeq ($(VCF_TYPE), multi-sample)
call_variants: $(MULTI_VCF)
else
call_variants: $(VCF).tbi
endif

# Run all steps from alignment to variant calling.
all: call_variants
> @echo "All steps completed."




