set -uex

#
# Keep MIT as a single copy
#

ASM=17XNL_SL_v3.fa
samtools faidx $ASM PY17XNL_MIT >mt.fa

samtools faidx mt.fa PY17XNL_MIT:2705-8659 >seq.fa
