#Instructions for using UGAP

-Dependencies that need to be installed by the user:

1. BLAST+ (must have blastn in PATH)
2. blastall (used by an external tool, but must be in PATH)
3. GATK (Recommend version < v3), must be in UGAP/bin directory, with UGAP being
   where UGAP is installed
4. SPAdes, should work with all versions, including 3.5
5. Musket [OPTIONAL] (only if you want to use this corrector)
6. genomeCoverageBed (part of BEDTOOLS)
7. bwa (must have version that supports BWA-MEM algorithm)
8. Samtools

-To run:

*You need a directory with paired-end reads in *.fastq.gz format. Names need to be similar
to what comes off Illumina sequencer (e.g. *R1_001.fastq.gz)

*You need to edit the "ugap_prep.py" script to reflect your installation location:

Change this line to fit your directory:

UGAP_PATH="/Users/jasonsahl/tools/UGAP"

*First you will need to generate a UGAP input file. You can do this like:

python ugap_prep.py -d reads -b /scratch/blastdb/nt -p 4 > ugap.config

*If you don't give a path to a BLAST database, it will ignore this functionality.

*Inspect the config file to make sure that you see commands that will be fed to UGAP

*Submit to the queue using the "all" partition in this case. This will change depending
on your own system

python ugap_univeral_controller.py -c ugap.config -o slurm -q all