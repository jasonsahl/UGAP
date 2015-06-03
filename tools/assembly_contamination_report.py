#!/usr/bin/env/python

"""Decon assemblies"""

import sys
import os
from optparse import OptionParser
import subprocess
from subprocess import Popen
try:
    from Bio import SeqIO
except:
    print "BioPython either not found or not configured correctly"
    sys.exit()
try:
    from ugap.util import *
except:
    print "Environment not set correctly, correct ugap_single.py environment"
    sys.exit()

def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print '%s file cannot be opened' % option
        sys.exit()

def bwa_single(genome_index,read_1,sam_file, processors, log_file='',**my_opts):
    mem_arguments = ['bwa', 'mem', '-v', '2', '-M', '-t', '%s' % processors]
    for opt in my_opts.items():
        mem_arguments.extend(opt)
    mem_arguments.extend([genome_index,read_1])

    if log_file:
       try:
           log_fh = open(log_file, 'w')
       except:
           print log_file, 'could not open'
    else:
        log_fh = PIPE

    try:
        sam_fh = open(sam_file, 'w')
    except:
        print sam_file, 'could not open'
    bwa = Popen(mem_arguments, stderr=log_fh, stdout=sam_fh)
    bwa.wait()

def bwa_paired(genome_index,read_1,read_2,sam_file, processors, log_file='',**my_opts):
    mem_arguments = ['bwa', 'mem', '-v', '2', '-M', '-t', '%s' % processors]
    for opt in my_opts.items():
        mem_arguments.extend(opt)
    mem_arguments.extend([genome_index,read_1,read_2])

    if log_file:
       try:
           log_fh = open(log_file, 'w')
       except:
           print log_file, 'could not open'
    else:
        log_fh = PIPE

    try:
        sam_fh = open(sam_file, 'w')
    except:
        print sam_file, 'could not open'
    bwa = Popen(mem_arguments, stderr=log_fh, stdout=sam_fh)
    bwa.wait()

def run_bwa(ref, read_1, read_2, name):
    read_group = '@RG\tID:%s\tSM:vac6wt\tPL:ILLUMINA\tPU:vac6wt' % name
    print "aligning reads to reference"
    if "NULL" in read_2:
        bwa_single(ref,read_1,"out.sam", processors, log_file='sam.log',**{'-R':read_group})
        print "alignment done"
    else:
        bwa_paired(ref, read_1, read_2, "out.sam", processors, log_file='sam.og',**{'-R':read_group})
        print "alignment done"
            
def main(forward_read,name,reverse_read,assembly,blast_nt):
    dependencies = ['bwa','samtools','genomeCoverageBed','blastn']
    for dependency in dependencies:
        ra = subprocess.check_call('which %s > /dev/null 2>&1' % dependency, shell=True)
        if ra == 0:
            pass
        else:
            print "%s is not in your path, but needs to be!" % dependency
            sys.exit()
    start_dir = os.getcwd()
    start_path = os.path.abspath("%s" % start_dir)
    forward_path = os.path.abspath("%s" % forward_read)
    if "NULL" not in reverse_read:
        reverse_path = os.path.abspath("%s" % reverse_read)
    blast_nt_path = os.path.abspath("%s" % blast_nt)
    assembly_path = os.path.abspath("%s" % assembly)
    if os.path.isfile("%s_renamed.bam" % name):
        print "will use pre-built bam file"
    else:
        print "building bam file"
        os.system("bwa index %s > /dev/null 2>&1" % assembly_path)
        if "NULL" not in reverse_read:
            run_bwa("%s" % assembly_path, forward_path, reverse_path, name)
        else:
            run_bwa("%s" % assembly_path, forward_path, "NULL", name)
    make_bam("%s.sam" % name, name)
    #All functions in the util.py file
    get_seq_length("%s" % assembly_path, name)
    subprocess.check_call("tr ' ' '\t' < %s.tmp.txt > %s.genome_size.txt" % (name, name), shell=True)
    get_coverage("%s_renamed.bam" % name,"%s.genome_size.txt" % name, name)
    remove_column("%s.tmp.out" % name, name)
    sum_coverage("%s.coverage.out" % name, 3, name)
    merge_files_by_column(0,"%s.genome_size.txt" % name, "%s.amount_covered.txt" % name, "%s.results.txt" % name)
    report_stats("%s.results.txt" % name, "%s_renamed_header.bam" % name, name)
    doc("%s.coverage.out" % name, "%s.genome_size.txt" % name, name, 3)
    os.system("cp %s_3_depth.txt %s/UGAP_assembly_results" % (name,start_path))
    sum_totals("%s_3_depth.txt" % name, name, "%s/UGAP_assembly_results/%s_coverage.txt" % (start_path,name))
    slice_assembly("%s" % assembly_path,200,"%s.chunks.fasta" % name)
    lengths = get_contig_lengths("%s" % assembly_path)
    print ""
    print "beginning BLAST"
    print ""
    subprocess.check_call("blastn -query %s.chunks.fasta -db %s -outfmt '7 std stitle' -dust no -evalue 0.01 -num_threads 4 -out %s.blast.out" % (name,blast_nt,name), shell=True)
    print "BLAST finished"
    print ""
    os.system("cp %s.blast.out %s_blast_report.txt" % (name,name))
    os.system("sort -u -k 1,1 %s.blast.out > %s.blast.uniques" % (name,name))
    merge_blast_with_coverages("%s.blast.uniques" % name, "%s_3_depth.txt" % name, lengths,name)
    os.system("sed 's/ /_/g' %s.depth_blast_merged.txt > %s.tmp.txt" % (name,name))
    os.system("sort -u -k 1,1 %s.tmp.txt | sort -gr -k 3,3 > %s_blast_depth_merged.txt" % (name,name))
    find_missing_coverages("%s_3_depth.txt" % name, "%s_blast_depth_merged.txt" % name, lengths,name)
    os.system("sort -u -k 1,1 %s.new.txt | sort -gr -k 5,5 > %s_blast_depth_merged.txt" % (name,name))
    #print "cleaning up"
    #os.system("rm sam.log tmp.txt genome_size.txt tmp.out coverage.out amount_covered.txt results.txt error.log blast.out new.txt")

if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-f", "--forward", dest="forward_read",
                      help="forward read, must be *.fastq.gz [REQUIRED]",
                      action="callback", callback=test_file, type="string")
    parser.add_option("-n", "--name", dest="name",
                      help="sample name [REQUIRED]",
                      action="store", type="string")
    parser.add_option("-v", "--reverse", dest="reverse_read",
                      help="reverse read, must be *.fastq.gz, CAN be missing",
                      action="store", type="string", default="NULL")
    parser.add_option("-a", "--assembly", dest="assembly",
                      help="genome assembly to decon [REQUIRED]",
                      action="callback", callback=test_file, type="string")
    parser.add_option("-b", "--blast_nt", dest="blast_nt",
                      help="path to genbank nt database [REQUIRED]",
                      action="store", type="string")

    options, args = parser.parse_args()
    mandatories = ["forward_read","name","assembly","blast_nt"]
    for m in mandatories:
        if not options.__dict__[m]:
            print "\nMust provide %s.\n" %m
            parser.print_help()
            exit(-1)
    main(options.forward_read,options.name,options.reverse_read,options.assembly,
         options.blast_nt)
