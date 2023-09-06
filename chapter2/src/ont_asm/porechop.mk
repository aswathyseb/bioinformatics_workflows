# Requirements and installation
# git clone https://github.com/rrwick/Porechop.git
# cd Porechop
# python3 setup.py install

#
# These commands trim any adapter sequences present in the raw nanopore reads.
#

# Nanopore reads
READS ?= reads/nano_dna.fq.gz

# Output folder
OUT_DIR ?= trimmed

# Trimmed output
TRIM_OUT = trimed/trimmed.fq

# Output folder
OUT_DIR = $(shell dirname ${TRIM_OUT})

#mkdir -p $OUT

# Makefile customizations.
.RECIPEPREFIX = >
.DELETE_ON_ERROR:
.ONESHELL:
MAKEFLAGS += --warn-undefined-variables --no-print-directory

# Print usage information.
usage::
> @echo "#"
> @echo "# porechop.mk: Remove adapter sequence using porechop"
> @echo "#"
> @echo "# make trim READS=${READS} TRIM_OUT=${TRIM_OUT}"
> @echo "#"

# Trim adapter sequences from Nanopore data
${TRIM_OUT}:${READS}
>mkdir -p ${OUT_DIR}
>porechop -i ${READS} -o ${TRIM_OUT} 1>${OUT_DIR}/porechop.log
#Generate stats before and after trimming
>seqkit stats ${READS} ${TRIM_OUT} >${OUT_DIR}/read_stats.txt

trim:${TRIM_OUT}
>ls -l ${TRIM_OUT}
