#!/usr/bin/env python

from __future__ import division
import random
import os
import logging
from Bio import SeqIO
import subprocess
from subprocess import Popen
import glob
from igs.threading import functional as p_func
import threading
import decimal
import re
import sys


def subsample_reads(input_fastq,output_fastq):
    #adapted from pythonforbiologists.com example
    from gzip import GzipFile
    import gzip
    import random
    """I changed this on 8/15/2016: Testing to see if this is enough or whether we should still increase"""
    number_to_sample = 20000000
    with GzipFile(input_fastq) as input:
        num_lines = sum([1 for line in input])
        total_records = int(num_lines / 4)
    if int(total_records)<int(number_to_sample):
        os.system("cp %s %s" % (input_fastq,output_fastq))
    else:
        outfile = gzip.open(output_fastq, "wb")
        output_sequence_sets = (set(random.sample(xrange(total_records + 1), number_to_sample)))
        record_number = 0
        with GzipFile(input_fastq) as input:
            for line1 in input:
                line2 = input.next()
                line3 = input.next()
                line4 = input.next()
                if record_number in output_sequence_sets:
                    outfile.write(line1)
                    outfile.write(line2)
                    outfile.write(line3)
                    outfile.write(line4)
                else:
                    pass
                record_number += 1
        outfile.close()

def report_stats(results, bam, name):
    infile = open(results, "rU")
    outfile = open("%s_breadth.txt" % name, "w")
    print >> outfile, name,"\n",
    for line in infile:
        fields = line.split()
        chromosome = fields[0]
        try:
            amount = (int(fields[2])/int(fields[1]))*100
            print >> outfile,chromosome,"\t",amount,"\n",
        except:
            print >> outfile, chromosome,"\t","0","\n",
            sys.exc_clear()
    infile.close()
    outfile.close()

def merge_files_by_column(column, file_1, file_2, out_file):
    """Takes 2 file and merge their columns based on the column. It is assumed
    that line ordering in the files do not match, so we read both files into memory
    and join them"""
    join_map = {}
    for line in open(file_1):
        line.strip()
        row = line.split()
        column_value = row.pop(column)
        join_map[column_value] = row
    for line in open(file_2):
        line.strip()
        row = line.split()
        column_value = row.pop(column)
        if column_value in join_map:
            join_map[column_value].extend(row)
    fout = open(out_file, 'w')
    for k, v in join_map.iteritems():
        fout.write('\t'.join([k] + v) + '\n')
    fout.close()

def doc(coverage, genome_size, name, suffix):
    incov = open(coverage, "U")
    ingenom = open(genome_size, "U")
    outfile = open("%s_%s_depth.txt" % (name, suffix), "w")
    all = [ ]
    my_dict = {}
    for line in incov:
        fields=line.split()
        fields = map(lambda s: s.strip(), fields)
        all.append(fields)
    for x, y in all:
        if int(y)>int(1):
           try:
                my_dict[x].append(y)
           except KeyError:
                my_dict[x] = [y]
        else:
           continue
    new_dict={}
    for k,v in my_dict.iteritems():
        ints = map(int, v)
        new_dict.update({k:sum(ints)})
    genome_size_dict = {}
    for line in ingenom:
        fields = line.split()
        genome_size_dict.update({fields[0]:fields[1]})
    print >> outfile, name,"\n",
    for k,v in new_dict.iteritems():
        print >> outfile, k,"\t",round(int(v)/int(genome_size_dict.get(k)),0)
    for y,z in genome_size_dict.iteritems():
        if y not in new_dict:
                print >> outfile, y,"\t","0"

def sum_coverage(coverage,cov,name):
    infile = open(coverage, "rU")
    outfile = open("%s.amount_covered.txt" % name, "w")
    all = [ ]
    dict = {}
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
    for k,v in dict.iteritems():
        print >> outfile, k+"\t"+str(len(v))
    infile.close()
    outfile.close()

def get_coverage(bam, size, name):
    """does the actual work"""
    subprocess.check_call("genomeCoverageBed -d -ibam %s -g %s > %s.tmp.out" % (bam,size,name), shell=True)

def remove_column(temp_file, name):
    infile = open(temp_file, "rU")
    outfile = open("%s.coverage.out" % name, "w")
    my_fields = [ ]
    for line in infile:
        fields=line.split()
        del fields[1]
        my_fields.append(fields)
    for x in my_fields:
        print >> outfile, "\t".join(x)
    infile.close()
    outfile.close()

def get_seq_length(ref, name):
    """uses BioPython in order to calculated the length of
    each fasta entry in the reference fasta"""
    infile = open(ref, "rU")
    outfile = open("%s.tmp.txt" % name, "w")
    for record in SeqIO.parse(infile, "fasta"):
        print >> outfile,record.id,len(record.seq)
    infile.close()
    outfile.close()

def run_trimmomatic(trim_path, processors, forward_path, reverse_path, ID, ugap_path, length):
    args=['java','-jar','%s' % trim_path, 'PE', '-threads', '%s' % processors,
              '%s' % forward_path, '%s' % reverse_path, '%s.F.paired.fastq.gz' % ID, 'F.unpaired.fastq.gz',
	      '%s.R.paired.fastq.gz' % ID, 'R.unpaired.fastq.gz', 'ILLUMINACLIP:%s/bin/illumina_adapters_all.fasta:4:30:10:1:true' % ugap_path,
	      'MINLEN:%s' % length]
    vcf_fh = open('%s.trimmomatic.out' % ID, 'w')
    log_fh = open('%s.trimmomatic.log' % ID, 'w')
    try:
        trim = Popen(args, stderr=vcf_fh, stdout=log_fh)
        trim.wait()
    except:
        log_isg.logPrint("problem encountered with trimmomatic")

def slice_assembly(infile, keep_length, outfile):
    """Keep length will be 200"""
    input=open(infile, "rU")
    output = open(outfile, "w")
    start=0
    #end=keep_length
    for record in SeqIO.parse(input,"fasta"):
        seqlength = len(record.seq)
        print >> output,">"+record.id+"\n",
        if seqlength<250:
            print >> output, record.seq[start:keep_length]
        elif seqlength<350:
            print >> output, record.seq[100:300]
        elif seqlength<450:
            print >> output, record.seq[200:400]
        else:
            print >> output, record.seq[200:400]

    input.close()
    output.close()

def find_missing_coverages(depth, merged, lengths, name):
    all_ids = {}
    outfile = open("%s.new.txt" % name, "w")
    for line in open(depth, "U"):
        fields = line.split()
        if len(fields)==1:
            pass
        else:
            all_ids.update({fields[0]:fields[1]})
    for k,v in all_ids.iteritems():
        hits = []
        nohits = []
        for line in open(merged, "U"):
            fields = line.split()
            if k == fields[0]:
                print >> outfile, line,
                hits.append("1")
            else:
                nohits.append("1")
        allhits = hits + nohits
        if len(nohits)==len(allhits):
            print >> outfile, str(k)+"\t"+"N/A"+"\t"+"no_blast_hit"+"\t"+str(lengths.get(k))+"\t"+str(v)
    outfile.close()

def merge_blast_with_coverages(blast_report, coverages, lengths, name):
    from operator import itemgetter
    coverage_dict = {}
    out_list = []
    outfile = open("%s.depth_blast_merged.txt" % name, "w")
    #dictionary contains contig:coverage
    for line in open(coverages, "U"):
        fields = line.split()
        if len(fields)==1:
            pass
        else:
            coverage_dict.update({fields[0]:fields[1]})
    for line in open(blast_report, "U"):
        file_list = []
        newline = line.strip()
        #out_list = []
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
        print >> outfile, "\t".join(alist)

def get_contig_lengths(in_fasta):
    length_dict = {}
    for record in SeqIO.parse(open(in_fasta, "U"), "fasta"):
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
                    print "#Could not determine forward/reverse read status for input file"
        else:
            baseName, read  = m.groups()[0], m.groups()[3]
            if read == "_R1":
                forward_reads[baseName] = infile
            elif read == "_R2":
                reverse_reads[baseName] = infile
            else:
                print "#Could not determine forward/reverse read status for input file "
                fileSets[file_name_before_ext] = infile
                num_single_readsets += 1
    for sample in forward_reads:
        if sample in reverse_reads:
            fileSets[sample] = [forward_reads[sample],reverse_reads[sample]] # store pair
            num_paired_readsets += 1
        else:
            fileSets[sample] = [forward_reads[sample]] # no reverse found
            num_single_readsets += 1
            logging.info('#Warning, could not find pair for read:' + forward_reads[sample])
    for sample in reverse_reads:
        if sample not in fileSets:
            fileSets[sample] = reverse_reads[sample] # no forward found
            num_single_readsets += 1
            logging.info('#Warning, could not find pair for read:' + reverse_reads[sample])

    if num_paired_readsets > 0:
        logging.info('Total paired readsets found:' + str(num_paired_readsets))
    if num_single_readsets > 0:
        logging.info('Total single reads found:' + str(num_single_readsets))

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
    for record in SeqIO.parse(open(fasta_in, "U"), "fasta"):
        seqrecords.append(record)
    output_handle=open(fasta_out, "w")
    SeqIO.write(seqrecords, output_handle, "fasta")
    output_handle.close()

def filter_seqs(fasta_in, keep, name):
    kept_sequences=[]
    for record in SeqIO.parse(open(fasta_in, "U"), "fasta"):
        if len(record.seq) >= int(keep):
            kept_sequences.append(record)
    output_handle = open("%s.%s.spades.assembly.fasta" % (name,keep), "w")
    SeqIO.write(kept_sequences, output_handle, "fasta")
    output_handle.close()

def bwa(reference,read_1,read_2,sam_file, processors, log_file='',**my_opts):
    mem_arguments = ['bwa', 'mem', '-v', '2', '-M', '-t', '%s' % processors]
    for opt in my_opts.items():
        mem_arguments.extend(opt)
    if "null" in read_2:
        mem_arguments.extend([reference,read_1])
    else:
        mem_arguments.extend([reference,read_1,read_2])
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

def run_bwa(read_1, read_2, processors, name, reference):
    read_group = '@RG\tID:%s\tSM:vac6wt\tPL:ILLUMINA\tPU:vac6wt' % name
    bwa(reference,read_1,read_2,"%s.sam" % name,processors,log_file='sam.log',**{'-R':read_group})

def make_bam(in_sam, name):
    subprocess.check_call("samtools view -h -b -S %s > %s.1.bam 2> /dev/null" % (in_sam, name), shell=True)
    subprocess.check_call("samtools view -u -h -F4 -o %s.2.bam %s.1.bam > /dev/null 2>&1" % (name,name), shell=True)
    subprocess.check_call("samtools view -h -b -q1 -F4 -o %s.3.bam %s.2.bam > /dev/null 2>&1" % (name,name), shell=True)
    subprocess.check_call("samtools sort %s.3.bam %s_renamed > /dev/null 2>&1" % (name,name), shell=True)
    subprocess.check_call("samtools index %s_renamed.bam > /dev/null 2>&1" % name, shell=True)
    subprocess.check_call("rm %s.1.bam %s.2.bam %s.3.bam" % (name,name,name), shell=True)

def run_gatk(reference, processors, name, gatk):
    args = ['java', '-jar', '%s' % gatk, '-T', 'UnifiedGenotyper',
            '-R', '%s' % reference, '-nt', '%s' % processors,'-S', 'silent',
            '-mbq', '17', '-ploidy', '1', '-out_mode', 'EMIT_VARIANTS_ONLY',
            '-stand_call_conf', '100', '-stand_emit_conf', '100', '-I', '%s_renamed.bam' % name,
            '-rf', 'BadCigar']
    try:
        vcf_fh = open('%s.gatk.out' % name, 'w')
    except:
        print 'could not open gatk file'
    try:
        log_fh = open('%s.gatk.log' % name, 'w')
    except:
        print 'could not open log file'
    gatk_run = Popen(args, stderr=log_fh, stdout=vcf_fh)
    gatk_run.wait()

def parse_vcf(vcf, coverage, proportion):
    vcf_in = open(vcf, "U")
    vcf_out = open("filtered.vcf", "w")
    to_fix = ()
    for line in vcf_in:
        if line.startswith('#'):
           continue
        fields=line.split()
        if "INFO" not in fields[0]:
            snp_fields=fields[9].split(':')
            if int(len(snp_fields))>2:
                prop_fields=snp_fields[1].split(',')
                if int(snp_fields[2])>=coverage:
                    if int(prop_fields[1])/int(snp_fields[2])>=float(proportion):
                        print >> vcf_out, line,
                        to_fix=((fields[0],fields[1],fields[4]),)+to_fix
            else:
                pass
        else:
            continue
    vcf_in.close()
    vcf_out.close()
    return to_fix

def rename_multifasta(fasta_in, prefix, fasta_out):
    """rename mutli-fasta to something meaningful"""
    rec=1
    handle = open(fasta_out, "w")
    for record in SeqIO.parse(open(fasta_in), "fasta"):
        print >> handle, ">"+prefix+"_"+str(autoIncrement())
        print >> handle, record.seq
    handle.close()

def rename_for_prokka(name_list):
    """17 allows some room for the unique contig number"""
    reduced_name_list = name_list[:17]
    new_name = "".join(reduced_name_list)
    return new_name

def sum_totals(input, name, output):
    outfile = open(output, "w")
    coverages = []
    for line in open(input, "U"):
        fields = line.split()
        if len(fields)<2:
            pass
        else:
            coverages.append(float(fields[1]))
    print >> outfile, name, sum(coverages)/len(coverages)
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

def run_single_loop(forward_path,reverse_path,name,error_corrector,processors,keep,start_path,reduce,careful,UGAP_PATH,TRIM_PATH,PICARD_PATH,PILON_PATH,blast_nt,cov_cutoff,phiX_filter):
    """This should, in theory, get rid of phiX if it is present, assuming it's not in the reference"""
    if "NULL" not in reduce:
        #Reads will be depleted in relation to a given reference
        rv = subprocess.call(['which', 'bam2fastq'])
        if rv == 0:
            #I can make this function more efficient
            run_bwa(forward_path, reverse_path, processors, name, reduce)
            os.system("samtools view -bS %s.sam > %s.bam 2> /dev/null" % (name,name))
            os.system("bam2fastq -o %s#.fastq --no-aligned %s.bam > %s.reduce_results.txt" % (name,name,name))
            os.system("gzip %s_1.fastq %s_2.fastq" % (name,name))
            subsample_reads("%s_1.fastq.gz" % name, "%s.F.tmp.fastq.gz" % name)
            subsample_reads("%s_2.fastq.gz" % name, "%s.R.tmp.fastq.gz" % name)
        else:
            print "to deplete reads, you need to have bam2fastq installed. Reads will not be depleted"
    """The Ks parameter is obsolete in Spade 3.7+, this will likely be removed in the future"""
    if int(get_sequence_length_dev(forward_path))<=200 and int(get_sequence_length_dev(forward_path))>=100:
        #Uses default K values, based on SPADes recs
        ks = "21,33,55,77"
    elif int(get_sequence_length_dev(forward_path))>200:
        ks = "21,33,55,77,99,127"
    elif int(get_sequence_length_dev(forward_path))<100:
        ks = "21,33"
    else:
        pass
    #Gets the sequence length independently for each genomes
    length = (int(get_sequence_length_dev(forward_path)/2))
    #Sub-sample reads to 4 million in each direction: at some point this will become a tunable parameter
    #Checkpoint 1: subsample reads
    if os.path.isfile("%s.F.tmp.fastq.gz" % name):
        pass
    else:
        subsample_reads(forward_path, "%s.F.tmp.fastq.gz" % name)
        subsample_reads(reverse_path, "%s.R.tmp.fastq.gz" % name)
    #If trimmomatic has already been run, don't run again, trimmomatic requires PAIRED reads
    #Checkpoint 2: trimmomatic, usearch
    if os.path.isfile("%s.F.paired.fastq.gz" % name):
        pass
    else:
        run_trimmomatic(TRIM_PATH, processors, "%s.F.tmp.fastq.gz" % name, "%s.R.tmp.fastq.gz" % name, name, UGAP_PATH, length)
        if phiX_filter == "T":
            #try:
            #    print "Removing phiX from reads with USEARCH"
            try:
                subprocess.check_call("gunzip %s.F.paired.fastq.gz %s.R.paired.fastq.gz > /dev/null 2>&1" % (name,name), shell=True)
            except:
                pass
            #subprocess.check_call("usearch -filter_phix %s.F.paired.fastq -reverse %s.R.paired.fastq -output %s.F.tmp.fastq -output2 %s.R.tmp.fastq > /dev/null 2>&1" % (name,name,name,name), shell=True)
            #subprocess.check_call("mv %s.F.tmp.fastq %s.F.paired.fastq" % (name,name), shell=True)
            #subprocess.check_call("mv %s.R.tmp.fastq %s.R.paired.fastq" % (name,name), shell=True)
            #subprocess.check_call("gzip %s.F.paired.fastq %s.R.paired.fastq" % (name,name), shell=True)
            #def uclust_sort(usearch):
            #    """sort with Usearch. Updated to V6"""
            #    devnull = open("/dev/null", "w")
            #    cmd = ["%s" % usearch,
            #           "-sortbylength", "all_gene_seqs.out",
            #           "-output", "tmp_sorted.txt"]
            #    subprocess.call(cmd,stdout=devnull,stderr=devnull)
            #    devnull.close()
            #devnull = open("/dev/null", "w")
            #os.system('usearch -filter_phix %s.F.paired.fastq -reverse %s.R.paired.fastq -output >(gzip > %s.F.tmp.fastq.gz) -output2 >(gzip > %s.R.tmp.fastq.gz)' % (name,name,name,name))
            cmd = ["usearch","-filter_phix","%s.F.paired.fastq" % name,"-reverse","%s.R.paired.fastq" % name,"-output","%s.F.tmp.fastq" % name,
                  "-output2","%s.R.tmp.fastq" % name]
            #print cmd
            #subprocess.call(cmd,stdout=devnull,stderr=devnull)
            subprocess.call(cmd)
            os.system("mv %s.F.tmp.fastq %s.F.paired.fastq" % (name,name))
            os.system("mv %s.R.tmp.fastq %s.R.paired.fastq" % (name,name))
            os.system("pigz *.paired.fastq")
            #devnull.close()
            #except:
            #    print "usearch9 required for phiX filtering...exiting"
            #    sys.exit()
    #This next section runs spades according to the input parameters
    #Checkpoint 3: Spades assembly
    if os.path.isfile("%s.spades.assembly.fasta" % name):
        pass
    else:
        if error_corrector=="hammer":
            if careful == "T":
                subprocess.check_call("spades.py -o %s.spades -t %s -k %s --cov-cutoff %s --careful -1 %s.F.paired.fastq.gz -2 %s.R.paired.fastq.gz  > /dev/null 2>&1" % (name,processors,ks,cov_cutoff,name,name), shell=True)
            else:
                subprocess.check_call("spades.py -o %s.spades -t %s -k %s --cov-cutoff %s -1 %s.F.paired.fastq.gz -2 %s.R.paired.fastq.gz  > /dev/null 2>&1" % (name,processors,ks,cov_cutoff,name,name), shell=True)
        else:
            if careful == "T":
                subprocess.check_call("spades.py --only-assembler -o %s.spades -t %s -k %s --cov-cutoff %s --careful -1 %s.F.paired.fastq.gz -2 %s.R.paired.fastq.gz  > /dev/null 2>&1" % (name,processors,ks,cov_cutoff,name,name), shell=True)
            else:
                subprocess.check_call("spades.py --only-assembler -o %s.spades -t %s -k %s --cov-cutoff %s -1 %s.F.paired.fastq.gz -2 %s.R.paired.fastq.gz  > /dev/null 2>&1" % (name,processors,ks,cov_cutoff,name,name),
        os.system("cp %s.spades/contigs.fasta %s.spades.assembly.fasta" % (name,name))
    #filters contigs by a user-defined length threshold, defaults to 200nts
    filter_seqs("%s.spades.assembly.fasta" % name, keep, name)
    #This uses biopython to pretty up the sequences, but not sure it would affect downstream usability
    clean_fasta("%s.%s.spades.assembly.fasta" % (name,keep),"%s_cleaned.fasta" % name)
    #Cleans up the names for downstream apps
    rename_multifasta("%s_cleaned.fasta" % name, name, "%s_renamed.fasta" % name)
    #Here I align reads to this new assembly: 1st instance of read alignment
    subprocess.check_call("bwa index %s_renamed.fasta > /dev/null 2>&1" % name, shell=True)
    #Index renamed.fasta for calling variants
    os.system("samtools faidx %s_renamed.fasta 2> /dev/null" % name)
    #run_bwa("%s.F.paired.fastq.gz" % name, "%s.R.paired.fastq.gz" % name, processors, name,"%s_renamed.fasta" % name)
    if "NULL" not in reduce:
        run_bwa("%s_1.fastq.gz" % name, "%s_2.fastq.gz" % name, processors, name, "%s_renamed.fasta" % name)
        subprocess.check_call("bwa mem -R '@RG\tID:${%s}\tSM:vac6wt\tPL:ILLUMINA\tPU:vac6wt' -v 2 -M -t %s %s_1.fastq.gz %s_2.fastq.gz | samtools view -uS - | samtools sort -@ '%s - '%s_renamed.bam'" % (name,processors,name,name,processors,name), shell=True)
    else:
        #align depleted reads if the reduced option is selected. This section is currently being tested
        subprocess.check_call("bwa mem -R '@RG\tID:${%s}\tSM:vac6wt\tPL:ILLUMINA\tPU:vac6wt' -v 2 -M -t %s %s.F.paired.fastq.gz  %s.F.paired.fastq.gz  | samtools view -uS - | samtools sort -@ '%s - '%s_renamed.bam'" % (name,processors,name,name,processors,name), shell=True)
        #run_bwa(forward_path, reverse_path, processors, name, "%s_renamed.fasta" % name)
        #make_bam("%s.sam" % name, name)
    print "running Pilon"
    try:
        os.system("java -jar %s --threads %s --fix all,amb --genome %s_renamed.fasta --bam %s_renamed.bam --output %s_pilon > /dev/null 2>&1" % (PILON_PATH,processors,name,name,name))
    except:
        print "problem running Pilon. Exiting...."
        #instead of exiting here, I could just change the name and keep going
        sys.exit()
    rename_multifasta("%s_pilon.fasta" % name, name, "%s_final_assembly.fasta" % name)
    filter_seqs("%s_final_assembly.fasta" % name, keep, name)
    #filters again by minimum length, output is named %s.%s.spades.assembly.fasta
    try:
        subprocess.check_call("sed -i 's/\\x0//g' %s.%s.spades.assembly.fasta" % (name,keep), shell=True, stderr=open(os.devnull, "w"))
    except:
        print "problem fixing missing spaces"
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
            #os.system("cp %s.%s.spades.assembly.fasta %s.prokka.fasta" % (name,keep,small_name))
            rename_multifasta("%s.%s.spades.assembly.fasta" % (name,keep), small_name, "%s.prokka.fasta" % small_name)
            os.system("prokka --prefix %s --locustag %s --centre %s --compliant --mincontiglen %s --strain %s %s.prokka.fasta > /dev/null 2>&1" % (small_name,small_name,small_name,keep,small_name,small_name))
        subprocess.check_call("cp %s/*.* %s/UGAP_assembly_results" % (name,start_path), shell=True, stderr=open(os.devnull, "w"))
    except:
        print "Prokka was not run, so no annotation files will be included"
    #Copies these files to your output directory, whether or not the previous commands were successful
    os.system("cp %s.%s.spades.assembly.fasta %s/UGAP_assembly_results/%s_final_assembly.fasta" % (name,keep,start_path,name))
    #I need to check the number of contigs, then decide whether or not to run bwa again
    os.system("bwa index %s.%s.spades.assembly.fasta > /dev/null 2>&1" % (name,keep))
    if "NULL" not in reduce:
        run_bwa("%s_1.fastq.gz" % name, "%s_2.fastq.gz" % name, processors, name, "%s.%s.spades.assembly.fasta" % (name,keep))
    else:
        run_bwa(forward_path,reverse_path,processors,name,"%s.%s.spades.assembly.fasta" % (name,keep))
    make_bam("%s.sam" % name, name)
    #This is for the per contig coverage routine. This can likely be replaced
    get_seq_length("%s.%s.spades.assembly.fasta" % (name,keep), name)
    subprocess.check_call("tr ' ' '\t' < %s.tmp.txt > %s.genome_size.txt" % (name, name), shell=True)
    get_coverage("%s_renamed.bam" % name,"%s.genome_size.txt" % name, name)
    remove_column("%s.tmp.out" % name, name)
    sum_coverage("%s.coverage.out" % name, 3, name)
    merge_files_by_column(0,"%s.genome_size.txt" % name, "%s.amount_covered.txt" % name, "%s.results.txt" % name)
    report_stats("%s.results.txt" % name, "%s_renamed_header.bam" % name, name)
    doc("%s.coverage.out" % name, "%s.genome_size.txt" % name, name, 3)
    os.system("cp %s_3_depth.txt %s/UGAP_assembly_results" % (name,start_path))
    sum_totals("%s_3_depth.txt" % name, name, "%s/UGAP_assembly_results/%s_coverage.txt" % (start_path,name))
    #End of section that can likely be replaced
    if "NULL" not in blast_nt:
        slice_assembly("%s.%s.spades.assembly.fasta" % (name,keep),int(keep),"%s.chunks.fasta" % name)
        lengths = get_contig_lengths("%s.%s.spades.assembly.fasta" % (name,keep))
        subprocess.check_call("blastn -task blastn -query %s.chunks.fasta -db %s -outfmt '7 std stitle' -dust no -evalue 0.01 -num_threads %s -out %s.blast.out" % (name, blast_nt, processors, name), shell=True)
        os.system("cp %s.blast.out %s/UGAP_assembly_results/%s_blast_report.txt" % (name, start_path, name))
        os.system("sort -u -k 1,1 %s.blast.out > %s.blast.uniques" % (name, name))
        merge_blast_with_coverages("%s.blast.uniques" % name, "%s_3_depth.txt" % name, lengths, name)
        os.system("sed 's/ /_/g' %s.depth_blast_merged.txt > %s.tmp.txt" % (name,name))
        os.system("sort -u -k 1,1 %s.tmp.txt | sort -gr -k 3,3 > %s/UGAP_assembly_results/%s_blast_depth_merged.txt" % (name, start_path, name))
        find_missing_coverages("%s_3_depth.txt" % name, "%s/UGAP_assembly_results/%s_blast_depth_merged.txt" % (start_path, name), lengths, name)
        os.system("sort -u -k 1,1 %s.new.txt | sort -gr -k 5,5 > %s/UGAP_assembly_results/%s_blast_depth_merged.txt" % (name, start_path, name))
    else:
        print "BLAST not run"
