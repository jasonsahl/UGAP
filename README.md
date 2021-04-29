### Instructions for using UGAP  

#### Dependencies for UGAP include:  

1. SPAdes  
2. bwa (must have version that supports BWA-MEM algorithm).  
3. Samtools  
4. BioPython  
5. blast+ (optional)  
6. bbmap  
7. seqtk  
8. bedtools  

### Installing with conda

```conda create -n ugap python=3.7```  
```conda activate ugap```  
```conda install -c conda-forge mamba```  
```mamba install -c bioconda spades bwa samtools biopython blast bbmap seqtk bedtools```  
```mamba install -c anaconda pigz```  
```git clone https://github.com/jasonsahl/UGAP.git```  

#### To run:

*You need a directory with paired-end reads in *.fastq.gz format. Names need to be similar to what comes off Illumina sequencer (e.g. *R1_001.fastq.gz). Single end support has not been thoroughly tested  

*You need to edit the "ugap_prep.py" script to reflect your installation location:  

Change this line to fit your directory:  

UGAP_PATH="/Users/jasonsahl/tools/UGAP"  

*First you will need to generate a UGAP input file. You can do this like:  

```python ugap_prep.py -d reads -b /scratch/blastdb/nt -p 4 > ugap.config```  

*If you don't give a path to a BLAST database, it will ignore this functionality and run the "sendsketch" method in bbmap instead  

*Inspect the config file to make sure that you see commands that will be fed to UGAP  

*Submit to the queue using the "all" partition in this case. This will change depending
on your own system  

```python ugap_univeral_controller.py -c ugap.config -o slurm -q all```  

