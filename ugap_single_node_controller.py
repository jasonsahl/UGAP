#!/usr/bin/env python


"""one controller to rule them all"""

import sys
import os
from optparse import OptionParser
from popen2 import popen2
import errno

try:
    from igs.utils import functional as func
    from igs.utils import logging
    from igs.threading import functional as p_func
except:
    print "Your environment is not set correctly.  Please correct your UGAP PATH and try again"
    sys.exit()


def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print '%s file cannot be opened' % option
        sys.exit()

def test_options(option, opt_str, value, parser):
    if "slurm" in value:
        setattr(parser.values, option.dest, value)
    elif "torque" in value:
        setattr(parser.values, option.dest, value)
    elif "sge" in value:
        setattr(parser.values, option.dest, value)
    else:
        print "option not supported, choose from slurm, torque, or sge"
        sys.exit()

def parse_config_file(config_file):
    datasets = ()
    infile = open(config_file, "U")
    for line in infile:
        if line.startswith("#"):
            pass
        else:
            fields = line.split()
            datasets=((fields[0],fields[1],fields[2],fields[3],fields[4],fields[5],fields[6],fields[7],fields[8],fields[9],fields[10],fields[11]),)+datasets
            processors = fields[7]
            ugap_path = fields[9]
    return datasets, processors, ugap_path

def run_pipeline_mixed(forward_path,reverse_path,name,error_corrector,processors,keep,coverage,proportion,start_path,reduce,careful, UGAP_PATH, TRIM_PATH, PICARD_PATH, PILON_PATH, GATK_PATH, blast_nt, cov_cutoff):
    if "NULL" not in reduce:
        #Reads will be depleted in relation to a given reference
        rv = subprocess.call(['which', 'bam2fastq'])
        if rx == 0:
            run_bwa("%s" % forward_path, "%s" % reverse_path, processors, name,"%s" % reduce)
            os.system("samtools view -bS %s.sam > %s.bam 2> /dev/null" % (name,name))
            os.system("bam2fastq -o %s#.fastq --no-aligned %s.bam > reduce_results.txt" % (name,name))
            os.system("gzip %s_1.fastq %s_2.fastq" % (name,name))
            os.system("cp %s_1.fastq.gz %s" % (name,forward_path))
            os.system("cp %s_2.fastq.gz %s" % (name,reverse_path))
        else:
            print "to deplete reads, you need to have bam2fastq installed. Reads will not be depleted"
    if int(get_sequence_length(forward_path, name))<=200 and int(get_sequence_length(forward_path, name))>=100:
        #Uses default K values, based on SPADes recs
        ks = "21,33,55,77"
    elif int(get_sequence_length(forward_path, name))>200:
        ks = "21,33,55,77,99,127"
    elif int(get_sequence_length(forward_path, name))<100:
        ks = "21,33"
    else:
        pass
    #Gets the sequence length independently for each genomes
    length = (int(get_sequence_length(forward_path,name)/2))
    #If trimmomatic has already been run, don't run again, trimmomatic requires PAIRED reads
    if os.path.isfile("%s.F.paired.fastq.gz" % name):
        pass
    else:
        run_trimmomatic(TRIM_PATH, processors, forward_path, reverse_path, name, UGAP_PATH, length)
    #This next section runs spades according to the input parameters
    if error_corrector=="hammer":
        if careful == "T":
            subprocess.check_call("spades.py -o %s.spades -t %s -k %s --cov-cutoff %s --careful -1 %s.F.paired.fastq.gz -2 %s.R.paired.fastq.gz  > /dev/null 2>&1" % (name,processors,ks,cov_cutoff,name,name), shell=True)
        else:
            subprocess.check_call("spades.py -o %s.spades -t %s -k %s --cov-cutoff %s -1 %s.F.paired.fastq.gz -2 %s.R.paired.fastq.gz  > /dev/null 2>&1" % (name,processors,ks,cov_cutoff,name,name), shell=True)
    elif error_corrector=="musket":
        ab = subprocess.call(['which', 'musket'])
        if ab == 0:
            pass
        else:
            print "musket isn't in your path, but needs to be!"
            sys.exit()
        #Need to test again with Musket
        subprocess.check_call("musket -k 17 8000000 -p %s -omulti %s -inorder %s.F.paired.fastq.gz %s.R.paired.fastq.gz > /dev/null 2>&1" % (processors,name,name,name), shell=True)
        subprocess.check_call("mv %s.0 %s.0.musket.fastq.gz" % (name,name), shell=True)
        subprocess.check_call("mv %s.1 %s.1.musket.fastq.gz" % (name,name), shell=True)
        if careful == "T":
            subprocess.check_call("spades.py -o %s.spades -t %s -k %s --cov-cutoff %s --only-assembler --careful -1  %s.0.musket.fastq.gz -2 %s.1.musket.fastq.gz > /dev/null 2>&1" % (name,processors,ks,cov_cutoff,name,name), shell=True)
        else:
            subprocess.check_call("spades.py -o %s.spades -t %s -k %s --cov-cutoff %s --only-assembler -1  %s.0.musket.fastq.gz -2 %s.1.musket.fastq.gz > /dev/null 2>&1" % (name,processors,ks,cov_cutoff,name,name), shell=True)
    else:
        if careful == "T":
            subprocess.check_call("spades.py -o %s.spades -t %s -k %s --cov-cutoff %s --only-assembler --careful -1 %s.F.paired.fastq.gz -2 %s.R.paired.fastq.gz > /dev/null 2>&1" % (name,processors,ks,cov_cutoff,name,name), shell=True)
        else:
            subprocess.check_call("spades.py -o %s.spades -t %s -k %s --cov-cutoff %s --only-assembler -1 %s.F.paired.fastq.gz -2 %s.R.paired.fastq.gz > /dev/null 2>&1" % (name,processors,ks,cov_cutoff,name,name), shell=True)
    os.system("cp %s.spades/contigs.fasta %s.spades.assembly.fasta" % (name,name))
    #filters contigs by a user-defined length threshold, defaults to 200nts
    filter_seqs("%s.spades.assembly.fasta" % name, keep, name)
    #This uses biopython to pretty up the sequences, but not sure it would affect downstream usability
    clean_fasta("%s.%s.spades.assembly.fasta" % (name,keep),"%s_cleaned.fasta" % name)
    #Cleans up the names for downstream apps
    rename_multifasta("%s_cleaned.fasta" % name, name, "%s_renamed.fasta" % name)
    """This section is actively being removed"""
    #Here I align reads to this new assembly
    #subprocess.check_call("bwa index %s_renamed.fasta > /dev/null 2>&1" % name, shell=True)
    #Index renamed.fasta for calling variants
    #os.system("samtools faidx %s_renamed.fasta 2> /dev/null" % name)
    #run_bwa("%s.F.paired.fastq.gz" % name, "%s.R.paired.fastq.gz" % name, processors, name,"%s_renamed.fasta" % name)
    #make_bam("%s.sam" % name, name)
    #os.system("java -jar %s/CreateSequenceDictionary.jar R=%s_renamed.fasta O=%s_renamed.dict > /dev/null 2>&1" % (PICARD_PATH, name, name))
    #run_gatk("%s_renamed.fasta" % name, processors, name, "%s" % GATK_PATH)
    """run_bam_coverage stuff here"""
    #os.system("java -jar %s/AddOrReplaceReadGroups.jar INPUT=%s_renamed.bam OUTPUT=%s_renamed_header.bam SORT_ORDER=coordinate RGID=%s RGLB=%s RGPL=illumina RGSM=%s RGPU=name CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT > /dev/null 2>&1" % (PICARD_PATH,name,name,name,name,name))
    #os.system("echo %s_renamed_header.bam > %s.bam.list" % (name,name))
    #This is redundant with per contig coverages discussed below
    #os.system("java -jar %s -R %s_renamed.fasta -T DepthOfCoverage -o %s_coverage -I %s.bam.list -rf BadCigar > /dev/null 2>&1" % (GATK_PATH,name,name,name))
    #os.system("samtools index %s_renamed_header.bam 2> /dev/null" % name)
    #process_coverage(name)
    #end of potentially redundant code
    #This next routine tries to fix SNPs if they are identified. Could be redundant with Pilon
    #try:
    #    to_fix=parse_vcf("%s.gatk.out" % name, coverage, proportion)
    #    log_isg.logPrint("number of SNPs to fix in %s = %s" % (name,len(to_fix)))
    #    if int(len(to_fix))>=1:
    #        try:
    #            fasta_to_tab("%s_renamed.fasta" % name, name)
    #            fix_assembly("%s.out.tab" % name, to_fix, name)
    #            os.system("cp %s_corrected_assembly.fasta %s_renamed.fasta" % (name,name))
    #        except:
    #            print "error correction failed for some reason"
    #    else:
    #        pass
    #except:
    #    print "couldn't correct SNPs, likely due to a GATK error"
    #    pass
    #Runs Pilon, I assume that it runs correctly
    os.system("java -jar %s --threads %s --genome %s_renamed.fasta --frags %s_renamed_header.bam --output %s_pilon > /dev/null 2>&1" % (PILON_PATH,processors,name,name,name))
    rename_multifasta("%s_pilon.fasta" % name, name, "%s_final_assembly.fasta" % name)
    filter_seqs("%s_final_assembly.fasta" % name, keep, name)
    #filters again by minimum length, output is named %s.%s.spades.assembly.fasta
    try:
        subprocess.check_call("sed -i 's/\\x0//g' %s.%s.spades.assembly.fasta" % (name,keep), shell=True, stderr=open(os.devnull, "w"))
    except:
        print "problem fixing missing spaces"
        pass
    clean_fasta("%s.%s.spades.assembly.fasta" % (name,keep),"%s/UGAP_assembly_results/%s_final_assembly.fasta" % (start_path,name))
    try:
        #Runs Prokka, if it's installed
        os.system("prokka --prefix %s --locustag %s --compliant --mincontiglen %s --strain %s %s.%s.spades.assembly.fasta > /dev/null 2>&1" % (name,name,keep,name,name,keep))
        subprocess.check_call("cp %s/*.* %s/UGAP_assembly_results" % (name,start_path), shell=True, stderr=open(os.devnull, "w"))
    except:
        print "Prokka was not run, so no annotation files will be included"
        pass
    #Copies these files to your output directory, whether or not the previous commands were successful
    os.system("cp %s.%s.spades.assembly.fasta %s/UGAP_assembly_results/%s_final_assembly.fasta" % (name,keep,start_path,name))
    #os.system("cp coverage_out.txt %s/UGAP_assembly_results/%s_coverage.txt" % (start_path,name))
    #I have to re-run bwa on the new assembly, as the contig lengths can change from the original SPAdes assembly
    os.system("bwa index %s.%s.spades.assembly.fasta > /dev/null 2>&1" % (name,keep))
    run_bwa("%s.F.paired.fastq.gz" % name, "%s.R.paired.fastq.gz" % name, processors, name,"%s.%s.spades.assembly.fasta" % (name,keep))
    make_bam("%s.sam" % name, name)
    #This is for the per contig coverage routine
    get_seq_length("%s.%s.spades.assembly.fasta" % (name,keep))
    subprocess.check_call("tr ' ' '\t' < tmp.txt > genome_size.txt", shell=True)
    get_coverage("%s_renamed.bam" % name,"genome_size.txt")
    remove_column("tmp.out")
    sum_coverage("coverage.out",coverage)
    merge_files_by_column(0,"genome_size.txt", "amount_covered.txt", "results.txt")
    report_stats("results.txt", "%s_renamed_header.bam" % name, name)
    doc("coverage.out", "genome_size.txt", name, coverage)
    os.system("cp %s_%s_depth.txt %s/UGAP_assembly_results" % (name,coverage,start_path))
    sum_totals("%s_%s_depth.txt" % (name,coverage), name, "%s/UGAP_assembly_results/%s_coverage.txt" % (start_path,name))
    if "NULL" not in blast_nt:
        slice_assembly("%s.%s.spades.assembly.fasta" % (name,keep),keep,"%s.chunks.fasta" % name)
        lengths = get_contig_lengths("%s.%s.spades.assembly.fasta" % (name,keep))
        print lengths
        subprocess.check_call("blastn -query %s.chunks.fasta -db %s -outfmt '7 std stitle' -dust no -evalue 0.01 -num_threads %s -out blast.out" % (name, blast_nt, processors), shell=True) 
        os.system("cp blast.out %s/UGAP_assembly_results/%s_blast_report.txt" % (start_path, name))
        os.system("sort -u -k 1,1 blast.out > blast.uniques")
        merge_blast_with_coverages("blast.uniques", "%s_%s_depth.txt" % (name,coverage), lengths)
        os.system("sed 's/ /_/g' depth_blast_merged.txt > tmp.txt")
        os.system("sort -u -k 1,1 tmp.txt | sort -gr -k 3,3 > %s/UGAP_assembly_results/%s_blast_depth_merged.txt" % (start_path, name))
        find_missing_coverages("%s_%s_depth.txt" % (name,coverage), "%s/UGAP_assembly_results/%s_blast_depth_merged.txt" % (start_path, name), lengths)
        os.system("sort -u -k 1,1 new.txt | sort -gr -k 5,5 > %s/UGAP_assembly_results/%s_blast_depth_merged.txt" % (start_path, name))
    else:
        print "BLAST not run"


def main(config_file, memory):
    datasets, processors, ugap_path=parse_config_file(config_file)
    UGAP_PATH=ugap_path
    PICARD_PATH=UGAP_PATH+"/bin/"
    TRIM_PATH=UGAP_PATH+"/bin/trimmomatic-0.30.jar"
    PILON_PATH=UGAP_PATH+"/bin/pilon-1.12.jar"
    start_dir = os.getcwd()
    start_path = os.path.abspath("%s" % start_dir)
    try:
        os.makedirs('%s/UGAP_assembly_results' % start_path)
        os.makedirs('%s/work_directory' % start_path)
    except OSError, e:
        if e.errno != errno.EEXIST:raise 
    if "NULL" != reduce:
        reduce_path=os.path.abspath("%s" % reduce)
    effective_jobs = int(int(memory)/8000)
    if effective_jobs <=1:
        effective_jobs = 1
    effective_processors = int(int(processors)/effective_jobs)
    os.chdir("%s/work_directory" % start_dir) 
    def _perform_workflow(data):
        tn, f = data
        run_single_loop(f[1],f[2],f[0],f[3],f[7],f[5],start_path,f[6],f[8],UGAP_PATH,TRIM_PATH,PICARD_PATH,PILON_PATH,f[10],f[11])
    results = set(p_func.pmap(_perform_workflow,
                              datasets,
                              num_workers=effective_jobs))

if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-c", "--config", dest="config_file",
                      help="config file that populates the UGAP single assembly",
                      action="callback", callback=test_file, type="string")
    parser.add_option("-m", "--memory", dest="memory",
                      help="amount of memory on the server, defaults to 48G, enter 48000",
                      action="store", type="string", default="48000")
    options, args = parser.parse_args()
    mandatories = ["config_file"]
    for m in mandatories:
        if not options.__dict__[m]:
            print "\nMust provide %s.\n" %m
            parser.print_help()
            exit(-1)
    main(options.config_file,options.memory)
        
