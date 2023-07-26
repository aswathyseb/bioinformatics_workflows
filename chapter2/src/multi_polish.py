import sys
#
# usage: 
# python multi_polish.py out_dir asm.fa R1 R2 REF >multi_polish.sh 

dirname = sys.argv[1]
REF = sys.argv[2]
R1 = sys.argv[3]
R2 = sys.argv[4]
GENOME = sys.argv[5]

CPU = 20

USE_FREEBAYES=True


# Read simulation parameters.
ERR = 0
N = 1000000
LEN = 150
READ1 = 'read1.fq'
READ2 = 'read2.fq'

# Nanopolish output
POLISH1 = 'step_2.fa'
# Homopolymer corrected output
POLISH2 = 'step_3.fa'
# Final correction.
POLISH3 = 'step_4.fa'

ALL1 = f"{POLISH1}.all.vcf.gz"
FILT1 = f"{POLISH1}.fb.sim.filt.vcf.gz"
# FILT1 = f"{POLISH1}.filt.vcf.gz"


setup = f"""
# Initialization

set -uex

rm -rf {dirname}
mkdir -p {dirname}
mkdir -p {dirname}/index
mkdir -p {dirname}/nextpolish
mkdir -p {dirname}/homopolymer_correction
cd {dirname}
cat ../{REF} > step_1.fa
echo > log.txt
"""

print(setup)

# Polish1 : Polish with nextPolish.
cmd_nextPolish = f"""
# *** Polish assmebly with nextPolish
cd nextpolish
ls ../../{R1} ../../{R2} >sgs.fofn
genome=../step_1.fa
echo -e "task = best\nworkdir=./rundir\ngenome = $genome\nsgs_fofn = sgs.fofn" > run.cfg
nextPolish run.cfg
cat ./rundir/genome.nextpolish.fasta >../{POLISH1}
samtools faidx ../{POLISH1}
cd ../
        """

print(cmd_nextPolish)

index1 = f"index/{POLISH1}"
bam_sim = f"{POLISH1}.sim.dna.bam"

# Polish2 : Homopolymer correction
cmd_correct = f"""
# *** Correct homopolymer erros.

cd homopolymer_correction
# simulate error-free reads.
pywgsim -S 0 -r 0 -R 0 -N {N} -e {ERR} -1 {LEN} -2 {LEN} ../../{GENOME} {READ1} {READ2} > log.txt 2>&1

# Build index
bwa index -p ../{index1} ../{POLISH1} >> log.txt  2>&1

# Align the simulated reads to assembly
bwa mem -t {CPU} ../{index1} {READ1} {READ2} 2>> log.txt | samtools sort > {bam_sim} 2>> log.txt
samtools index {bam_sim}

# Call variants with freebayes.
freebayes -f ../{POLISH1} {bam_sim} -C 5 -p 1 | bcftools norm --threads {CPU} -f ../{POLISH1} -o {POLISH1}.fb.sim.vcf
bgzip {POLISH1}.fb.sim.vcf
bcftools index {POLISH1}.fb.sim.vcf.gz

# Filter long indels
cat {POLISH1}.fb.sim.vcf.gz |bcftools filter -Oz --IndelGap 5 -i 'TYPE!~"snp" && abs(strlen(ALT)-strlen(REF))>=5  \
&& FORMAT/AO >5 && GT=="1" ' >{FILT1}
bcftools index  {FILT1}

# Consensus call, creating new build
cat ../{POLISH1} | bcftools consensus -H 1 {FILT1} > {POLISH2}
cp {POLISH2} ../
#cp {POLISH1}.fb.vcf.gz* ../
#cp {FILT1}* ../
cd ../
"""
print(cmd_correct)

# Polish3 : Correction with illumina reads.

STEPS = 4

for x in range(3, STEPS + 1):

    ref = f"step_{x}.fa"
    bam_dna = f"{ref}.illumina.dna.bam"
    all = f"{ref}.vars.all.vcf.gz"
    vcf = f"{ref}.vars.trusted.vcf.gz"
    index = f"index/{ref}"
    cmd1 = f'''

#
#
# *** alignment step {x}
#

bwa index -p {index} {ref}  >> log.txt  2>&1

# Illumina
bwa mem -t {CPU} {index} ../{R1} ../{R2} 2>> log.txt | samtools sort > {bam_dna} 2>> log.txt

# Call variants with freebayes.
freebayes -f {ref} {bam_dna} -C 5 -p 1 |  bcftools norm --threads {CPU} -f {ref} -o {ref}.fb.vcf 
bgzip {ref}.fb.vcf

# Filter variants
cat {ref}.fb.vcf.gz |bcftools filter -Oz --IndelGap 5 -i  'AF!=0 && GT=="1"' >{ref}.fb.filt.vcf.gz
bcftools index {ref}.fb.filt.vcf.gz


# Index and cleanup
ls -1 *.bam | parallel samtools index {{}}
rm -f *.ht2 *.ann *.bwt *.amb
    '''
    print(cmd1)

    if x != STEPS:
        seqid = f"{x + 1}"
        new = f"step_{seqid}.fa"

        cmd_freebayes = f'''
        # *** consensus call, creating build {seqid}
        cat {ref} | bcftools consensus -H 1 {ref}.fb.filt.vcf.gz > {new}
        '''

        cmd_bcftools = f'''
    # *** consensus call, creating build {seqid}
    cat {ref} | bcftools consensus -H 1 {vcf} > {new}
            '''

        if USE_FREEBAYES:
            print(cmd_freebayes)
        else:
            print(cmd_bcftools)

final_cmd=f"""cp {new} assembly_polished.fa
    
cd ../"""
print(final_cmd)

