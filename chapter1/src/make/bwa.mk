# The reference genome.
REF ?= refs/genome.fa

# The sample name.
SAMPLE ?= S1

# Library type/Alignment mode.
LIB ?= PE

# Number of CPUs to use.
NCPU ?= 4

# BWA flags to use.
BFLAGS ?=

#Sam flags to use.
SFLAGS ?=

# The genome index folder
IDX_DIR ?= $(dir ${REF})/idx

# Genome index prefix
IDX = ${IDX_DIR}/$(notdir ${REF})

# A bwa index file.
IDX_FILE = ${IDX}.amb

# The input reads.
R1 ?= data/${SAMPLE}_R1.fq
R2 ?= data/${SAMPLE}_R2.fq

# The alignment file.
BAM ?= results/bam/${SAMPLE}.bam

# The alignment directory
ALIGN_DIR ?= $(dir ${BAM})

# Makefile customizations.
.RECIPEPREFIX = >
.DELETE_ON_ERROR:
.ONESHELL:
MAKEFLAGS += --warn-undefined-variables --no-print-directory

#Print usage.
bwa_usage:
> @echo "#"
> @echo "# Aligns data to the reference using bwa."
> @echo "#"
>@echo "make align REF=${REF}  LIB=${LIB} SAMPLE=${SAMPLE} R1=${R1} R2=${R2} BAM"
> @echo "#"
> @echo "# REF    : Reference Genome"
> @echo "# LIB    : Library type, PE or SE"
> @echo "# SAMPLE : Sample name"
> @echo "# R1     : Read1"
> @echo "# R2     : Read2"
> @echo "#"

# Rule to create bwa index.
${IDX_FILE}: ${REF}
>@mkdir -p ${IDX_DIR}
>bwa index -p ${IDX} ${REF}

# Align in paired end mode.
ifeq ($(LIB), PE)
${BAM}: ${R1} ${R2} ${IDX_FILE}
> mkdir -p ${ALIGN_DIR}
> bwa mem ${BFLAGS} -t ${NCPU} ${IDX} ${R1} ${R2} | samtools view -h ${SFLAGS}| samtools sort -@ ${NCPU} >${BAM}
endif

# Single end alignment.
ifeq ($(LIB), SE)
${BAM}: ${R1} ${IDX_FILE}
> mkdir -p ${ALIGN_DIR}
> bwa mem ${BFLAGS} -t ${NCPU} ${IDX} ${R1} | samtools view -h ${SFLAGS}  | samtools sort -@ ${NCPU} >${BAM}
endif

# Create the BAM index file.
${BAM}.bai: ${BAM}
> samtools index ${BAM}

# The index target depends on the IDX file.
index: ${IDX_FILE}
>ls -l ${IDX_FILE}

# The align target depends on the BAM file.
bwa_align: ${BAM}.bai
>ls -l ${BAM}

align: index bwa_align
>ls -l ${IDX_FILE}
>ls -l ${BAM}

