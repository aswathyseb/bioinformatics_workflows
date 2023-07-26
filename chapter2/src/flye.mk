# Requirements and installation
# conda install flye
#
# De novo assembly with Flye assembler.
#

# Nanopore reads
READS ?= reads/nano_dna.fq.gz

# Read type
# Options are  nano-corr,nano-raw,nano-hq 
RTYPE ?= nano-raw

# Output folder
OUT ?= flye_out

# Estimated genome size
GSIZE ?= 23m

# Assembled fasta
ASM = ${OUT}/assembly.fasta

# Number of CPUs to use
CPU = 12

# Other Flye specific flags
#FLAGS = --scaffold
FLAGS =

# Makefile customizations.
.RECIPEPREFIX = >
.DELETE_ON_ERROR:
.ONESHELL:
MAKEFLAGS += --warn-undefined-variables --no-print-directory

# Print usage information.
usage::
> @echo "#"
> @echo "# flye.mk: Create genome assembly using Flye asssembler"
> @echo "#"
> @echo "# make assembly READS=${READS} GENOME_SIZE=${GSIZE} OUT=${OUT} RTYPE=${RTYPE}"
> @echo "#"

#Flye assembler
#
${ASM}:${READS}
>mkdir -p ${OUT}
>flye --nano-corr ${READS} --out-dir ${OUT} --genome-size ${GSIZE} -t ${CPU} ${FLAGS}

assembly:${ASM}
>ls -l ${ASM}
