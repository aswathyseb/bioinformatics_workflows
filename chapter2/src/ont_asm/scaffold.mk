#
# These commands scaffold an assembly based on a reference using ragtag
#


# Draft assembly
ASM ?= assembly.fasta

# Reference genome
REF ?= reference.fasta

# Output folder
OUT ?= ragtag_output

# Number of CPUs to use
CPU = 12

# Other specific flags
FLAGS = 

# Scaffolded fasta file
SCAFF = ${OUT}/ragtag.scaffold.fasta

# Print usage information.
usage::
	@echo "#"
	@echo "# scaffold.mk: Scaffold an assembly based on a reference using ragtag.py"
	@echo "#"
	@echo "# make scaffold ASM=${ASM} REF=${REF} OUT=${OUT}"
	@echo "#"

# Scaffold
${SCAFF}: ${READS} ${REF}
	ragtag.py scaffold ${REF} ${ASM} -o ${OUT}

scaffold:${SCAFF}
	ls -l ${SCAFF}
