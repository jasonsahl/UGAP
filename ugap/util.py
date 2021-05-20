#!/usr/bin/env python

from __future__ import division
import random
import os
import logging
from Bio import SeqIO
import subprocess
from subprocess import Popen
import glob
import threading
import decimal
import re
import sys

def subsample_reads_dev(input_fastq,output_fastq):
    #Trying out seqtk
    from gzip import GzipFile
    import gzip
    number_to_sample = 4000000
    with GzipFile(input_fastq) as input:
        num_lines = sum([1 for line in input])
        total_records = int(num_lines / 4)
        if (int(total_records)-5000)<int(number_to_sample):
            os.system("cp %s %s" % (input_fastq,output_fastq))
        else:
            subprocess.check_call("seqtk sample -s100 %s 4000000 | gzip > tmp.fastq.gz" % input_fastq,stdout=open(os.devnull, 'wb'),stderr=open(os.devnull, 'wb'), shell=True)
            os.system("mv tmp.fastq.gz %s" % output_fastq)

def report_stats(results, bam, name):
    outfile = open("%s_breadth.txt" % name, "w")
    outfile.write(name)
    outfile.write("\n")
    with open(results) as infile:
        for line in infile:
            fields = line.split()
            chromosome = fields[0]
            try:
                amount = (int(fields[2])/int(fields[1]))*100
                outfile.write(chromosome+"\t"+str(amount)+"\n")
            except:
                outfile.write(str(chromosome)+"\t"+"0"+"\n")
    outfile.close()

def merge_files_by_column(column, file_1, file_2, out_file):
    """Takes 2 file and merge their columns based on the column. It is assumed
    that line ordering in the files do not match, so we read both files into memory
    and join them"""
    join_map = {}
    with open(file_1) as my_file1:
        for line in my_file1:
            row = line.split()
            column_value = row.pop(column)
            join_map[column_value] = row
    with open(file_2) as my_file2:
        for line in open(file_2):
            row = line.split()
            column_value = row.pop(column)
            if column_value in join_map:
                join_map[column_value].extend(row)
    fout = open(out_file,'w')
    for k, v in join_map.items():
        fout.write('\t'.join([k] + v) + '\n')
    fout.close()

def doc(coverage, genome_size, name, suffix):
    outfile = open("%s_depth.txt" % name, "w")
    my_dict = {}
    with open(coverage) as my_cov:
        for line in my_cov:
            fields= line.split()
            fields = list(map(lambda s: s.strip(), fields))
            if int(fields[1])>1:
                try:
                    my_dict[fields[0]].append(int(fields[1]))
                except KeyError:
                    my_dict[fields[0]] = [int(fields[1])]
    genome_size_dict = {}
    with open(genome_size) as ingenom:
        for line in ingenom:
            fields = line.split()
            genome_size_dict.update({fields[0]:fields[1]})
    outfile.write(name+"\n")
    for k,v in my_dict.items():
        outfile.write(str(k)+"\t"+str(round(int(sum(v))/int(genome_size_dict.get(k))))+"\n")
    for y,z in genome_size_dict.items():
        if y not in my_dict:
            outfile.write(str(y)+"\t"+"0"+"\t"+"\n")
    outfile.close()

def sum_coverage(coverage,cov,name):
    outfile = open("%s.amount_covered.txt" % name, "w")
    all = []
    dict = {}
    with open(coverage) as infile:
        for line in infile:
            fields=line.split()
            fields = map(lambda s: s.strip(), fields)
            all.append(fields)
    for x, y in all:
        if int(y)>int(cov):
            try:
               dict[x].append(y)
            except KeyError:
               dict[x] = [y]
        else:
            pass
    for k,v in dict.items():
        outfile.write(str(k)+"\t"+str(len(v))+"\n")
    outfile.close()

def get_coverage_dev(bam, size, name):
    """trying out a variation, reducing the number of dependencies"""
    subprocess.check_call("samtools depth -aa %s > %s.tmp.out" % (bam,name), stdout=open(os.devnull, 'wb'),stderr=open(os.devnull, 'wb'),shell=True)

def remove_column(temp_file, name):
    outfile = open("%s.coverage.out" % name, "w")
    my_fields = []
    with open(temp_file) as infile:
        for line in infile:
            fields=line.split()
            del fields[1]
            my_fields.append(fields)
    for x in my_fields:
        outfile.write("\t".join(x)+"\n")
    outfile.close()

def get_seq_length(ref, name):
    """uses BioPython in order to calculated the length of
    each fasta entry in the reference fasta"""
    outfile = open("%s.tmp.txt" % name, "w")
    with open(ref) as infile:
        for record in SeqIO.parse(infile, "fasta"):
            outfile.write(str(record.id)+"\t"+str(len(record.seq))+"\n")
    outfile.close()

def slice_assembly(infile, keep_length, outfile):
    output = open(outfile, "w")
    start=0
    with open(infile) as input:
        for record in SeqIO.parse(input,"fasta"):
            seqlength = len(record.seq)
            output.write(">"+str(record.id)+"\n")
            if seqlength<250:
                output.write(str(record.seq[start:200])+"\n")
            elif seqlength<350:
                output.write(str(record.seq[100:300])+"\n")
            elif seqlength<450:
                output.write(str(record.seq[200:400])+"\n")
            else:
                output.write(str(record.seq[200:400])+"\n")
    output.close()

def find_missing_coverages(depth, merged, lengths, name):
    all_ids = {}
    outfile = open("%s.new.txt" % name,"w")
    with open(depth) as my_file:
        for line in my_file:
            fields = line.split()
            if len(fields)==1:
                pass
            else:
                all_ids.update({fields[0]:fields[1]})
    for k,v in all_ids.items():
        hits = []
        nohits = []
        with open(merged) as my_file:
            for line in my_file:
                fields = line.split()
                if k == fields[0]:
                    outfile.write(line)
                    hits.append("1")
                else:
                    nohits.append("1")
        allhits = hits + nohits
        if len(nohits)==len(allhits):
            outfile.write(str(k)+"\t"+"N/A"+"\t"+"no_blast_hit"+"\t"+str(lengths.get(k))+"\t"+str(v)+"\n")
    outfile.close()

def merge_blast_with_coverages(blast_report,coverages,lengths,name):
    from operator import itemgetter
    coverage_dict = {}
    out_list = []
    outfile = open("%s.depth_blast_merged.txt" % name,"w")
    #dictionary contains contig:coverage
    with open(coverages) as my_file:
        for line in my_file:
            fields = line.split()
            if len(fields)==1:
                pass
            else:
                coverage_dict.update({fields[0]:fields[1]})
    with open(blast_report) as my_blast:
        for line in my_blast:
            file_list = []
            newline = line.strip()
            if line.startswith("#"):
                pass
            else:
                fields = newline.split("\t")
                single_list = []
                single_list.append(str(fields[0]))
                single_list.append(str(fields[12]))
                single_list.append(str(fields[10]))
                single_list.append(str(lengths.get(fields[0])))
                single_list.append(str(coverage_dict.get(fields[0])))
                out_list.append(single_list)
    for alist in out_list:
        outfile.write("\t".join(alist))
        outfile.write("\n")
    outfile.close()

def merge_sendsketch_with_coverages(blast_report,coverages,lengths,name):
    from operator import itemgetter
    coverage_dict = {}
    out_list = []
    outfile = open("%s.depth_sendsketch_merged.txt" % name,"w")
    #dictionary contains contig:coverage
    with open(coverages) as my_file:
        for line in my_file:
            fields = line.split()
            if len(fields)==1:
                pass
            else:
                coverage_dict.update({fields[0]:fields[1]})
    with open(blast_report) as my_blast:
        for line in my_blast:
            file_list = []
            newline = line.strip()
            fields = newline.split("\t")
            out_list.append([fields[0],fields[1],fields[2],fields[3],str(coverage_dict.get(fields[0]))])
    for alist in out_list:
        outfile.write("\t".join(alist))
        outfile.write("\n")
    outfile.close()

def get_contig_lengths(in_fasta):
    length_dict = {}
    with open(in_fasta) as my_fasta:
        for record in SeqIO.parse(my_fasta,"fasta"):
           length_dict.update({record.id:len(record.seq)})
    return length_dict

def get_readFile_components(full_file_path):
    (file_path,file_name) = os.path.split(full_file_path)
    m1 = re.match("(.*).gz",file_name)
    ext = ""
    if m1 != None:
        ext = ".gz"
        file_name = m1.groups()[0]
    (file_name_before_ext,ext2) = os.path.splitext(file_name)
    full_ext = ext2+ext
    return(file_path,file_name_before_ext,full_ext)

#def read_file_sets(dir_path):
#    """match up pairs of sequence data, adapted from
#    https://github.com/katholt/srst2 will be tough to test
#    with variable names and read paths"""
#    fileSets = {}
#    forward_reads = {}
#    reverse_reads = {}
#    num_paired_readsets = 0
#    num_single_readsets = 0
#    for infile in glob.glob(os.path.join(dir_path, "*.fastq.gz")):
#        (file_path,file_name_before_ext,full_ext) = get_readFile_components(infile)
#        m=re.match("(_R*)", file_name_before_ext)
#        if m is None:
#            #m=re.match("(.*)("+"_R1"+")(_.*)$",file_name_before_ext)
#            m=re.match("(.*)("+"_R1"+")(.*)$",file_name_before_ext)
#            if m is not None:
#                (baseName,read) = m.groups()[0], m.groups()[1]
#                forward_reads[baseName] = infile
#            else:
                #m=re.match("(.*)("+"_R2"+")(_.*)$",file_name_before_ext)
#                m=re.match("(.*)("+"_R2"+")(.*)$",file_name_before_ext)
#                if m is not None:
#                    (baseName,read) = m.groups()[0], m.groups()[1]
#                    reverse_reads[baseName] = infile
#                else:
#                    fileSets[file_name_before_ext] = infile
#                    num_single_readsets += 1
                    #print("Could not determine forward/reverse read status for input file %s" % infile)
#        else:
#            baseName, read  = m.groups()[0], m.groups()[3]
#            if read == "_R1":
#                forward_reads[baseName] = infile
#            elif read == "_R2":
#                reverse_reads[baseName] = infile
#            else:
                #print("Could not determine forward/reverse read status for input file")
#                fileSets[file_name_before_ext] = infile
#                num_single_readsets += 1
#    for sample in forward_reads:
#        if sample in reverse_reads:
#            fileSets[sample] = [forward_reads[sample],reverse_reads[sample]] # store pair
#            num_paired_readsets += 1
##            fileSets[sample] = [forward_reads[sample]] # no reverse found
#            num_single_readsets += 1
            #logging.info('Warning, could not find pair for read:' + forward_reads[sample])
#            fileSets[sample] = reverse_reads[sample] # no forward found
#            num_single_readsets += 1
            #logging.info('Warning, could not find pair for read:' + reverse_reads[sample])
    #if num_paired_readsets > 0:
    #    logging.info('Total paired readsets found:' + str(num_paired_readsets))
    #if num_single_readsets > 0:
    #    logging.info('Total single reads found:' + str(num_single_readsets))
#    return fileSets

def read_file_sets(dir_path):
    fileSets = {}
    forward_reads = {}
    reverse_reads = {}
    num_paired_readsets = 0
    num_single_readsets = 0
    for infile in glob.glob(os.path.join(dir_path, "*.fastq.gz")):
        (file_path,file_name_before_ext,full_ext) = get_readFile_components(infile)
        m=re.match("(.*)(_S.*)(_L.*)(_R.*)(_.*)", file_name_before_ext)
        if m==None:
            m=re.match("(.*)("+"_R1"+")(_.*)$",file_name_before_ext)
            if m!=None:
                (baseName,read) = m.groups()[0], m.groups()[1]
                forward_reads[baseName] = infile
            else:
                m=re.match("(.*)("+"_R2"+")(_.*)$",file_name_before_ext)
                if m!=None:
                    (baseName,read) = m.groups()[0], m.groups()[1]
                    reverse_reads[baseName] = infile
                else:
                    pass
                    #print("#Could not determine forward/reverse read status for input file")
        else:
            baseName, read  = m.groups()[0], m.groups()[3]
            if read == "_R1":
                forward_reads[baseName] = infile
            elif read == "_R2":
                reverse_reads[baseName] = infile
            else:
                #print("#Could not determine forward/reverse read status for input file ")
                fileSets[file_name_before_ext] = infile
                num_single_readsets += 1
    for sample in forward_reads:
        if sample in reverse_reads:
            fileSets[sample] = [forward_reads[sample],reverse_reads[sample]] # store pair
            num_paired_readsets += 1
        else:
            fileSets[sample] = [forward_reads[sample]] # no reverse found
            num_single_readsets += 1
    for sample in reverse_reads:
        if sample not in fileSets:
            fileSets[sample] = reverse_reads[sample] # no forward found
            num_single_readsets += 1
    return fileSets

def get_seq_name(in_fasta):
    """used for renaming the sequences"""
    return os.path.basename(in_fasta)

def get_sequence_length_dev(fastq_in):
    from itertools import islice
    from gzip import GzipFile
    with GzipFile("%s" % fastq_in) as file:
        head = list(islice(file, 2))
    return len(head[1])

def clean_fasta(fasta_in, fasta_out):
    seqrecords=[]
    with open(fasta_in) as my_fasta:
        for record in SeqIO.parse(my_fasta, "fasta"):
            seqrecords.append(record)
    output_handle=open(fasta_out, "w")
    SeqIO.write(seqrecords, output_handle, "fasta")
    output_handle.close()

def filter_seqs(fasta_in, keep, name):
    kept_sequences=[]
    with open(fasta_in) as my_fasta:
        for record in SeqIO.parse(my_fasta, "fasta"):
            if len(record.seq) >= int(keep):
                kept_sequences.append(record)
    output_handle = open("%s.%s.spades.assembly.fasta" % (name,keep), "w")
    SeqIO.write(kept_sequences, output_handle, "fasta")
    output_handle.close()

def rename_multifasta(fasta_in, prefix, fasta_out):
    """rename mutli-fasta to something meaningful"""
    rec=1
    handle = open(fasta_out, "w")
    with open(fasta_in) as my_fasta:
        for record in SeqIO.parse(my_fasta, "fasta"):
            handle.write(">"+str(prefix)+"_"+str(autoIncrement())+"\n")
            handle.write(str(record.seq)+"\n")
    handle.close()

def rename_for_prokka(name_list):
    """17 allows some room for the unique contig number"""
    reduced_name_list = name_list[:17]
    new_name = "".join(reduced_name_list)
    return new_name

def sum_totals(input, name, output):
    outfile = open(output, "w")
    coverages = []
    with open(input) as my_file:
        for line in my_file:
            fields = line.split()
            if len(fields)<2:
                pass
            else:
                coverages.append(float(fields[1]))
    outfile.write(str(name)+"\t"+str(sum(coverages)/len(coverages))+"\n")
    outfile.close()

rec=1

def autoIncrement():
    global rec
    pStart = 1
    pInterval = 1
    if (rec == 0):
        rec = pStart
    else:
        rec += pInterval
        return rec

def run_sendsketch(fasta_in,name):
    from Bio import SeqIO
    import os
    import subprocess
    outfile = open("sendsketch_parsed.txt", "w")
    with open(fasta_in) as my_fasta:
        for record in SeqIO.parse(my_fasta, "fasta"):
            tmpfile = open("tmp.file.fasta", "w")
            tmpfile.write(">"+str(record.id)+"\n"+str(record.seq)+"\n")
            tmpfile.close()
            subprocess.check_call("sendsketch.sh overwrite=true address=nt in=tmp.file.fasta out=%s.sendsketch.out" % record.id,stdout=open(os.devnull, 'wb'),stderr=open(os.devnull, 'wb'), shell=True)
            with open("%s.sendsketch.out" % record.id) as my_file:
                try:
                    read_line = my_file.readlines()[3].strip()
                    fields = read_line.split("\t")
                    desired_fields = [record.id,fields[0],str(len(record.seq)),fields[10]]
                    outfile.write("\t".join(desired_fields)+"\n")
                except:
                    desired_fields = [record.id,"N/A",str(len(record.seq)),"no_hits"]
                    outfile.write("\t".join(desired_fields)+"\n")
            os.system("rm tmp.file.fasta %s.sendsketch.out" % record.id)
    outfile.close()

def run_single_loop(assembler,forward_path,reverse_path,name,error_corrector,processors,keep,start_path,reduce,careful,UGAP_PATH,
    TRIM_PATH,PICARD_PATH,PILON_PATH,blast_nt,cov_cutoff,phiX_filter,sample_type):
    if "NULL" not in reduce:
        #Reads will be depleted in relation to a given reference
        rv = subprocess.call(['which', 'bam2fastq'])
        if rv == 0:
            #I can make this function more efficient
            os.system("cp %s ./to_reduce.fasta" % reduce)
            os.system("bwa index to_reduce.fasta > /dev/null 2>&1")
            run_bwa("%s" % forward_path, "%s" % reverse_path, processors, name, "to_reduce.fasta")
            os.system("samtools index %s_renamed.bam" % name)
            os.system("bam2fastq -o %s#.fastq --no-aligned %s_renamed.bam > %s.reduce_results.txt" % (name,name,name))
            os.system("gzip %s_1.fastq %s_2.fastq" % (name,name))
            subsample_reads_dev("%s_1.fastq.gz" % name, "%s.F.tmp.fastq.gz" % name)
            subsample_reads_dev("%s_2.fastq.gz" % name, "%s.R.tmp.fastq.gz" % name)
        else:
            print("to deplete reads, you need to have bam2fastq installed. Reads will not be depleted")
    """The Ks parameter is obsolete in Spade 3.7+, this will likely be removed in the future"""
    #if int(get_sequence_length_dev(forward_path))<=200 and int(get_sequence_length_dev(forward_path))>=100:
        #Uses default K values, based on SPADes recs
        #ks = "21,33,55,77"
    #elif int(get_sequence_length_dev(forward_path))>200:
        #ks = "21,33,55,77,99,127"
    #elif int(get_sequence_length_dev(forward_path))<100:
        #ks = "21,33"
    #This is a placeholder for the K value required by spadees
    ks = "auto"
    #length is required for the trimming step with bbduk
    length = (int(get_sequence_length_dev(forward_path)/2))
    #Sub-sample reads to 4 million in each direction: at some point this will become a tunable parameter
    """Checkpoint 1: subsample reads"""
    if os.path.isfile("%s.F.tmp.fastq.gz" % name):
        pass
    else:
        if sample_type == "PE":
            subsample_reads_dev(forward_path, "%s.F.tmp.fastq.gz" % name)
            subsample_reads_dev(reverse_path, "%s.R.tmp.fastq.gz" % name)
        else:
            subsample_reads_dev(forward_path, "%s.F.tmp.fastq.gz" % name)
    #If trimmomatic has already been run, don't run again, trimmomatic requires PAIRED reads
    """Checkpoint 2: trimmomatic, usearch"""
    if os.path.isfile("%s.F.paired.fastq.gz" % name):
        pass
    else:
        if sample_type == "PE":
            #subprocess.check_call("bbduk.sh in=%s.F.tmp.fastq.gz in2=%s.R.tmp.fastq.gz ref=%s/bin/univec.fasta out=%s.F.paired.fastq.gz out2=%s.R.paired.fastq.gz minlen=%s overwrite=true" % (name,name,UGAP_PATH,name,name,length), shell=True)
            subprocess.check_call("bbduk.sh k=17 in=%s.F.tmp.fastq.gz in2=%s.R.tmp.fastq.gz ref=%s/bin/illumina_adapters_all.fasta out=%s.F.paired.fastq.gz out2=%s.R.paired.fastq.gz minlen=%s overwrite=true" % (name,name,UGAP_PATH,name,name,length), shell=True)
        elif sample_type == "SE":
            #subprocess.check_call("bbduk.sh in=%s.F.tmp.fastq.gz ref=%s/bin/univec.fasta out=%s.F.paired.fastq.gz minlen=%s overwrite=true" % (name,UGAP_PATH,name,length), shell=True)
            subprocess.check_call("bbduk.sh k=17 in=%s.F.tmp.fastq.gz ref=%s/bin/illumina_adapters_all.fasta out=%s.F.paired.fastq.gz minlen=%s overwrite=true" % (name,UGAP_PATH,name,length), shell=True)
        if phiX_filter == "T":
            try:
                if sample_type == "PE":
                    subprocess.check_call("gunzip %s.F.paired.fastq.gz %s.R.paired.fastq.gz > /dev/null 2>&1" % (name,name), shell=True)
                elif sample_type == "SE":
                    subprocess.check_call("gunzip %s.F.paired.fastq.gz > /dev/null 2>&1" % name, shell=True)
            except:
                pass
            if sample_type == "PE":
                cmd = ["usearch","-filter_phix","%s.F.paired.fastq" % name,"-reverse","%s.R.paired.fastq" % name,"-output","%s.F.tmp.fastq" % name,
                      "-output2","%s.R.tmp.fastq" % name]
            elif sample_type == "SE":
                cmd = ["usearch","-filter_phix","%s.F.paired.fastq" % name,"-output","%s.F.tmp.fastq" % name]
            try:
                devnull = open("/dev/null", "w")
                subprocess.call(cmd, stderr=devnull, stdout=devnull)
                try:
                    os.system("mv %s.F.tmp.fastq %s.F.paired.fastq" % (name,name))
                    os.system("mv %s.R.tmp.fastq %s.R.paired.fastq" % (name,name))
                    os.system("pigz *.paired.fastq")
                except:
                    os.system("mv %s.F.tmp.fastq %s.F.paired.fastq" % (name,name))
                    os.system("pigz *.paired.fastq")
            except:
                print("usearch9 required for phiX filtering...exiting")
                sys.exit()
    #This next section runs spades according to the input parameters
    #Checkpoint 3: Spades assembly
    if os.path.isfile("%s.spades.assembly.fasta" % name):
        pass
    elif os.path.isfile("%s.skesa.fasta" % name):
        pass
    else:
        if error_corrector=="hammer":
            if assembler=="spades":
                if careful == "T":
                    if sample_type == "PE":
                        subprocess.check_call("spades.py -o %s.spades -t %s -k %s --cov-cutoff %s --careful -1 %s.F.paired.fastq.gz -2 %s.R.paired.fastq.gz  > /dev/null 2>&1" % (name,processors,ks,cov_cutoff,name,name), shell=True)
                    elif sample_type == "SE":
                        subprocess.check_call("spades.py -o %s.spades -t %s -k %s --cov-cutoff %s --careful -s %s.F.paired.fastq.gz > /dev/null 2>&1" % (name,processors,ks,cov_cutoff,name), shell=True)
                else:
                    if sample_type == "PE":
                        subprocess.check_call("spades.py -o %s.spades -t %s -k %s --cov-cutoff %s -1 %s.F.paired.fastq.gz -2 %s.R.paired.fastq.gz  > /dev/null 2>&1" % (name,processors,ks,cov_cutoff,name,name), shell=True)
                    elif sample_type == "SE":
                        subprocess.check_call("spades.py -o %s.spades -t %s -k %s --cov-cutoff %s -s %s.F.paired.fastq.gz > /dev/null 2>&1" % (name,processors,ks,cov_cutoff,name), shell=True)
            #TODO: add SE support for skesa
            elif assembler=="skesa":
                subprocess.check_call("spades.py --only-error-correction -o %s.spades -t %s -1 %s.F.paired.fastq.gz -2 %s.R.paired.fastq.gz  > /dev/null 2>&1" % (name,processors,name,name), shell=True)
                """This workflow needs to be tested"""
                subprocess.check_call("skesa --gz --fastq %s.spades/corrected/%s.F.paired.fastq.00.0_0.cor.fastq.gz,%s.spades/corrected/%s.R.paired.fastq.00.0_0.cor.fastq.gz --cores %s --contigs_out %s.skesa.fasta > /dev/null 2>&1" % (name,name,name,name,processors,name), shell=True)
        else:
            if assembler=="spades":
                if careful == "T":
                    if sample_type == "PE":
                        subprocess.check_call("spades.py --only-assembler --careful -o %s.spades -t %s -k %s --cov-cutoff %s -1 %s.F.paired.fastq.gz -2 %s.R.paired.fastq.gz  > /dev/null 2>&1" % (name,processors,ks,cov_cutoff,name,name), shell=True)
                    elif sample_type == "SE":
                        subprocess.check_call("spades.py --only-assembler --careful -o %s.spades -t %s -k %s --cov-cutoff %s -s %s.F.paired.fastq.gz > /dev/null 2>&1" % (name,processors,ks,cov_cutoff,name), shell=True)
                else:
                    if sample_type == "PE":
                        subprocess.check_call("spades.py --only-assembler -o %s.spades -t %s -k %s --cov-cutoff %s -1 %s.F.paired.fastq.gz -2 %s.R.paired.fastq.gz  > /dev/null 2>&1" % (name,processors,ks,cov_cutoff,name,name), shell=True)
                    elif sample_type == "SE":
                        subprocess.check_call("spades.py --only-assembler -o %s.spades -t %s -k %s --cov-cutoff %s -s %s.F.paired.fastq.gz > /dev/null 2>&1" % (name,processors,ks,cov_cutoff,name), shell=True)
            else:
                #TODO: add SE support for skesa
                subprocess.check_call("skesa --gz --fastq %s.F.paired.fastq.gz,%s.R.paired.fastq.gz --cores %s --contigs_out %s.skesa.fasta > /dev/null 2>&1" % (name,name,processors,name), shell=True)
    """need to rename stuff here, copying over the Skesa files if they exist"""
    if assembler=="spades":
        os.system("cp %s.spades/contigs.fasta %s.spades.assembly.fasta" % (name,name))
        filter_seqs("%s.spades.assembly.fasta" % name, keep, name)
        #This uses biopython to pretty up the sequences, but not sure it would affect downstream usability
        clean_fasta("%s.%s.spades.assembly.fasta" % (name,keep),"%s_cleaned.fasta" % name)
    else:
        filter_seqs("%s.skesa.fasta" % name, keep, name)
        clean_fasta("%s.%s.spades.assembly.fasta" % (name,keep),"%s_cleaned.fasta" % name)
    #Cleans up the names for downstream apps
    rename_multifasta("%s_cleaned.fasta" % name, name, "%s_renamed.fasta" % name)
    #Here I align reads to this new assembly: 1st instance of read alignment
    subprocess.check_call("bwa index %s_renamed.fasta" % name,stdout=open(os.devnull, 'wb'),stderr=open(os.devnull, 'wb'),shell=True)
    #Index renamed.fasta for calling variants
    #os.system("samtools faidx %s_renamed.fasta 2> /dev/null" % name)
    try:
        if "NULL" not in reduce:
            if sample_type == "PE":
                subprocess.check_call("bwa mem -v 2 -M -t 4 %s_renamed.fasta %s_1.fastq.gz %s_2.fastq.gz | samtools sort -l 0 -@ 4 - | samtools view -Su -o %s_renamed.bam -" % (name,name,name,name),stdout=open(os.devnull, 'wb'),stderr=open(os.devnull, 'wb'),shell=True)
            elif sample_type == "SE":
                subprocess.check_call("bwa mem -v 2 -M -t 4 %s_renamed.fasta %s_1.fastq.gz | samtools sort -l 0 -@ 4 - | samtools view -Su -o %s_renamed.bam -" % (name,name,name),stdout=open(os.devnull, 'wb'),stderr=open(os.devnull, 'wb'),shell=True)
            os.system("samtools index %s_renamed.bam" % name)
        else:
            #align depleted reads if the reduced option is selected. This section is currently being tested
            #TODO: add support here
            if sample_type == "PE":
                #print("bwa mem -v 2 -M -t 4 %s_renamed.fasta %s %s | samtools sort -l 0 -@ 4 - | samtools view -Su -o %s_renamed.bam -" % (name,forward_path,reverse_path,name))
                subprocess.check_call("bwa mem -v 2 -M -t 4 %s_renamed.fasta %s %s | samtools sort -l 0 -@ 4 - | samtools view -Su -o %s_renamed.bam -" % (name,forward_path,reverse_path,name),stdout=open(os.devnull, 'wb'),stderr=open(os.devnull, 'wb'),shell=True)
            elif sample_type == "SE":
                subprocess.check_call("bwa mem -v 2 -M -t 4 %s_renamed.fasta %s | samtools sort -l 0 -@ 4 - | samtools view -Su -o %s_renamed.bam -" % (name,forward_path,name),stdout=open(os.devnull, 'wb'),stderr=open(os.devnull, 'wb'),shell=True)
            os.system("samtools index %s_renamed.bam" % name)
    except:
        print("problems running bwa")
        sys.exit()
    print("running Pilon")
    try:
        os.system("java -jar %s --threads %s --fix all,amb --genome %s_renamed.fasta --bam %s_renamed.bam --output %s_pilon > /dev/null 2>&1" % (PILON_PATH,processors,name,name,name))
    except:
        print("problem running Pilon. Exiting....")
        #instead of exiting here, I could just change the name and keep going
        sys.exit()
    rename_multifasta("%s_pilon.fasta" % name,name,"%s_final_assembly.fasta" % name)
    filter_seqs("%s_final_assembly.fasta" % name,keep,name)
    #filters again by minimum length, output is named %s.%s.spades.assembly.fasta
    try:
        subprocess.check_call("sed -i 's/\\x0//g' %s.%s.spades.assembly.fasta" % (name,keep), shell=True, stderr=open(os.devnull, "w"))
    except:
        print("problem fixing missing spaces")
    clean_fasta("%s.%s.spades.assembly.fasta" % (name,keep),"%s/UGAP_assembly_results/%s_final_assembly.fasta" % (start_path,name))
    try:
        #Runs Prokka, if it's installed
        name_chars = []
        for letter in name:
            name_chars.append(letter)
        if len(name_chars)<=20:
            os.system("prokka --prefix %s --locustag %s --centre %s --compliant --mincontiglen %s --strain %s %s.%s.spades.assembly.fasta > /dev/null 2>&1" % (name,name,name,keep,name,name,keep))
        else:
            """Tries to fix the short read limitation. Need to test"""
            small_name = rename_for_prokka(name_chars)
            rename_multifasta("%s.%s.spades.assembly.fasta" % (name,keep), small_name, "%s.prokka.fasta" % small_name)
            subprocess.check_call("prokka --prefix %s --locustag %s --centre %s --compliant --mincontiglen %s --strain %s %s.prokka.fasta" % (small_name,small_name,small_name,keep,small_name,small_name),
            stdout=open(os.devnull, 'wb'),stderr=open(os.devnull, 'wb'),shell=True)
        subprocess.check_call("cp %s/*.* %s/UGAP_assembly_results" % (name,start_path), shell=True, stderr=open(os.devnull, "w"))
    except:
        print("Prokka was not run or failed, so no annotation files will be included")
    #Copies these files to your output directory, whether or not the previous commands were successful
    os.system("cp %s.%s.spades.assembly.fasta %s/UGAP_assembly_results/%s_final_assembly.fasta" % (name,keep,start_path,name))
    #I need to check the number of contigs, then decide whether or not to run bwa again
    subprocess.check_call("bwa index %s.%s.spades.assembly.fasta" % (name,keep),stdout=open(os.devnull, 'wb'),stderr=open(os.devnull, 'wb'),shell=True)
    if "NULL" not in reduce:
        if sample_type == "PE":
            subprocess.check_call("bwa mem -v 2 -M -t 4 %s.%s.spades.assembly.fasta %s_1.fastq.gz %s_2.fastq.gz | samtools sort -l 0 -@ 4 - | samtools view -Su -o %s_renamed.bam -" % (name,keep,name,name,name),stdout=open(os.devnull, 'wb'),stderr=open(os.devnull, 'wb'),shell=True)
        elif sample_type == "SE":
            subprocess.check_call("bwa mem -v 2 -M -t 4 %s.%s.spades.assembly.fasta %s_1.fastq.gz | samtools sort -l 0 -@ 4 - | samtools view -Su -o %s_renamed.bam -" % (name,keep,name,name),stdout=open(os.devnull, 'wb'),stderr=open(os.devnull, 'wb'),shell=True)
    else:
        if sample_type == "PE":
            subprocess.check_call("bwa mem -v 2 -M -t 4 %s.%s.spades.assembly.fasta %s %s | samtools sort -l 0 -@ 4 - | samtools view -Su -o %s_renamed.bam -" % (name,keep,forward_path,reverse_path,name),stdout=open(os.devnull, 'wb'),stderr=open(os.devnull, 'wb'),shell=True)
        elif sample_type == "SE":
            subprocess.check_call("bwa mem -v 2 -M -t 4 %s.%s.spades.assembly.fasta %s | samtools sort -l 0 -@ 4 - | samtools view -Su -o %s_renamed.bam -" % (name,keep,forward_path,name),stdout=open(os.devnull, 'wb'),stderr=open(os.devnull, 'wb'),shell=True)
    #This is for the per contig coverage routine.
    get_seq_length("%s.%s.spades.assembly.fasta" % (name,keep), name)
    subprocess.check_call("tr ' ' '\t' < %s.tmp.txt > %s.genome_size.txt" % (name, name), shell=True)
    #get_coverage_dev("%s_renamed.bam" % name,"%s.genome_size.txt" % name, name)
    subprocess.check_call("samtools depth -aa %s_renamed.bam > %s.tmp.out" % (name,name), stdout=open(os.devnull, 'wb'),stderr=open(os.devnull, 'wb'),shell=True)
    remove_column("%s.tmp.out" % name, name)
    sum_coverage("%s.coverage.out" % name, 3, name)
    merge_files_by_column(0,"%s.genome_size.txt" % name,"%s.amount_covered.txt" % name,"%s.results.txt" % name)
    report_stats("%s.results.txt" % name, "%s_renamed_header.bam" % name, name)
    """The 3 suffix here is totally arbitrary and should be changed"""
    doc("%s.coverage.out" % name, "%s.genome_size.txt" % name, name, 3)
    os.system("cp %s_depth.txt %s/UGAP_assembly_results" % (name,start_path))
    sum_totals("%s_depth.txt" % name, name, "%s/UGAP_assembly_results/%s_coverage.txt" % (start_path,name))
    #End of section that can likely be replaced
    try:
        if "NULL" in blast_nt:
            """This means that sendsketch will be run"""
            lengths = get_contig_lengths("%s.%s.spades.assembly.fasta" % (name,keep))
            run_sendsketch("%s.%s.spades.assembly.fasta" % (name,keep),name)
            os.system("cp sendsketch_parsed.txt %s/UGAP_assembly_results/%s_sendsketch.txt" % (start_path,name))
            merge_sendsketch_with_coverages("sendsketch_parsed.txt","%s_depth.txt" % name,lengths,name)
            os.system("sed 's/ /_/g' %s.depth_sendsketch_merged.txt > %s.tmp.txt" % (name,name))
            os.system("sort -gr -k 5,5 %s.tmp.txt > %s/UGAP_assembly_results/%s.depth_sendsketch_merged.txt" % (name, start_path, name))
        else:
            slice_assembly("%s.%s.spades.assembly.fasta" % (name,keep),int(keep),"%s.chunks.fasta" % name)
            lengths = get_contig_lengths("%s.%s.spades.assembly.fasta" % (name,keep))
            subprocess.check_call("blastn -task blastn -query %s.chunks.fasta -db %s -outfmt '7 std stitle' -dust no -evalue 0.01 -num_threads %s -out %s.blast.out" % (name, blast_nt, processors, name), shell=True)
            os.system("cp %s.blast.out %s/UGAP_assembly_results/%s_blast_report.txt" % (name, start_path, name))
            os.system("sort -u -k 1,1 %s.blast.out > %s.blast.uniques" % (name, name))
            merge_blast_with_coverages("%s.blast.uniques" % name, "%s_depth.txt" % name, lengths, name)
            os.system("awk '{print $NF,$0}' %s.depth_blast_merged.txt | sort -nr | cut -f2- -d' ' > %s/UGAP_assembly_results/%s_blast_depth_merged.txt" % (name,start_path,name))
            find_missing_coverages("%s_depth.txt" % name, "%s/UGAP_assembly_results/%s_blast_depth_merged.txt" % (start_path, name), lengths, name)
    except:
        print("BLAST not run")
