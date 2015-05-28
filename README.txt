#Instructions for using UGAP (see Manual for more detailed instructions)

-Dependencies that MUST be installed by the user:

1. SPAdes, must be version 3.5
2. genomeCoverageBed (part of BEDTOOLS). Linux version included in /bin folder. If
   incompatible with your version, you must replace this version with your binary.
3. bwa (must have version that supports BWA-MEM algorithm)
4. Samtools
5. BioPython

-Optional dependencies for full functionality

1. BLAST+ (must have blastn in PATH). This is needed for troubleshooting mixtures
2. Musket (only if you want to correct reads with this tool)


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
