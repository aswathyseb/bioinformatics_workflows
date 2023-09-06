## Chapter2 - Assembly and annotation of _Plasmodium yoelli 17XNL_ strain

The worklows and scripts used in the assembly of P.yoelii 17XNL strain and annotation of its gene models are here.

The workflow Makefiles are

1. Makefile.ont - used in hybrid genome assembly
2. Makefile.pacbio - used in pacbio genome assembly
3. Makefile.ont - contains common operations like mapping, variant calling etc.

All the accessory scripts are in `src` folder.

**Nanopore-Illumina hybrid genome assembly**


    make -f Makefile.ont 
    #
    # Makefile.ont.mk: Create Nanopore-Illumina hybrid genome assembly
    #
    # make all READS=reads/nano_dna.fq.gz GENOME_SIZE=23m R1=reads/illumina_R1.fq.gz R2=reads/illumina_R2.fq.gz REF=refs/Py17X_Genome.fasta
    #
    # READS : Nanopore reads
    # R1,R2 : Illumina reads
    # GENOME_SIZE : Estimated genome size
    # REF: Reference genome (used for polishing)
    # 
    # The final output is assembly_polished.fa in the polish directory.

The folder `src/ont_asm`  contains the additional scripts that was used in making the final assembly.


**Pacbio genome assembly**

    make -f Makefile.pacbio 
    #
    # Makefile.pacbio.mk: Create Pacbio genome assembly with canu
    #
    # make all READS=reads GENOME_SIZE=23m REF=refs/Py17X_Genome.fasta
    #
    # READS : Data folder with pacbio subreads
    # GENOME_SIZE : Estimated genome size
    # REF: Reference genome
    #    
    # The final output is selected contigs based on reference alignment.
 

The folder `src/pacbio_asm`  contains the additional scripts that was used in making the final assembly.
   
    

