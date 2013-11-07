#!/usr/bin/env python

"""the ultimate genome assembly pipeline"""

from optparse import OptionParser
import sys
import os
import glob
import re
import logging
from igs.utils import logging as log_isg
from igs.threading import functional as p_func
import subprocess
from subprocess import Popen
from Bio import SeqIO
import errno
import threading

UGAP_PATH="/home/jsahl/UGAP"
os.system("export UGAP_DIR=%s" % UGAP_PATH)
os.system("export PYTHONPATH=$UGAP_DIR:$PYTHONPATH")
os.system("export PATH=$UGAP_DIR/bin:$PATH")
GATK_PATH=UGAP_PATH+"/bin/GenomeAnalysisTK.jar"
PICARD_PATH=UGAP_PATH+"/bin/CreateSequenceDictionary.jar"
TRIM_PATH=UGAP_PATH+"/bin/trimmomatic-0.30.jar"

def test_dir(option, opt_str, value, parser):
    if os.path.exists(value):
        setattr(parser.values, option.dest, value)
    else:
        print "directory of fastqs cannot be found"
        sys.exit()

def test_options(option, opt_str, value, parser):
    if "hammer" in value:
        setattr(parser.values, option.dest, value)
    elif "musket" in value:
        setattr(parser.values, option.dest, value)
    elif "none" in value:
        setattr(parser.values, option.dest, value)
    else:
        print "select from hammer, musket, or none"
        sys.exit()

def test_truths(option, opt_str, value, parser):
    if "T" in value:
        setattr(parser.values, option.dest, value)
    elif "F" in value:
        setattr(parser.values, option.dest, value)
    else:
        print "must select from T or F"
        sys.exit()

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
            m=re.match("(.*)("+"_1"+")$",file_name_before_ext)
            if m!=None:
                (baseName,read) = m.groups()
                forward_reads[baseName] = infile
            else:
                m=re.match("(.*)("+"_2"+")$",file_name_before_ext)
                if m!=None:
                    (baseName,read) = m.groups()
                    reverse_reads[baseName] = infile
                else:
                    print "Could not determine forward/reverse read status for input file"
        else:
            baseName, read  = m.groups()[0], m.groups()[3]
            if read == "_R1":
                forward_reads[baseName] = infile
            elif read == "_R2":
                reverse_reads[baseName] = infile
            else:
                print "Could not determine forward/reverse read status for input file "
                fileSets[file_name_before_ext] = infile
                num_single_readsets += 1
    for sample in forward_reads:
        if sample in reverse_reads:
            fileSets[sample] = [forward_reads[sample],reverse_reads[sample]] # store pair
            num_paired_readsets += 1
        else:
            fileSets[sample] = [forward_reads[sample]] # no reverse found
            num_single_readsets += 1
            logging.info('Warning, could not find pair for read:' + forward_reads[sample])
    for sample in reverse_reads:
        if sample not in fileSets:
            fileSets[sample] = reverse_reads[sample] # no forward found
            num_single_readsets += 1
            logging.info('Warning, could not find pair for read:' + reverse_reads[sample])
                                
    if num_paired_readsets > 0:
        logging.info('Total paired readsets found:' + str(num_paired_readsets))        
    if num_single_readsets > 0:
        logging.info('Total single reads found:' + str(num_single_readsets))

    return fileSets

def get_sequence_length(fastq_in, name):
    os.system("zcat %s > %s.tmp.fastq" % (fastq_in,name))
    lines = [ ]
    with open("%s.tmp.fastq" % name, "U") as f:
        for line in f.readlines()[1:2]:
            for x in line:
                lines.append(x)
    length=len(lines)
    os.system("rm -rf %s.tmp.fastq" % name)
    return length
    
                    
def run_image(assembly, name, length):
    my_k = int(length)-20
    args = ['image.pl', '-scaffolds','%s' % assembly,
            '-prefix', '%s' % name,'-iteration','1',
            '-all_iteration','2','-dir_prefix','%s' % name,
            '-kmer','%s' % my_k]
    try:
        vcf_fh = open('%s.image.out' % name, 'w')
    except:
        log_isg.logPrint('could not open image log file')
    try:
        log_fh = open('%s.image.log' % name, 'w')
    except:
        log_isg.logPrint('could not open log file')
    try:
        image = Popen(args, stderr=vcf_fh, stdout=log_fh)
        image.wait()
    except:
        log_isg.logPrint("problem encountered with image")

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
    subprocess.check_call("samtools view -h -b -S %s -o 1.bam 2> /dev/null" % in_sam, shell=True)
    subprocess.check_call("samtools view -u -h -F4 -o 2.bam 1.bam", shell=True)
    subprocess.check_call("samtools view -h -b -q1 -F4 -o 3.bam 2.bam", shell=True)
    subprocess.check_call("samtools sort 3.bam %s_pagit" % name, shell=True)
    subprocess.check_call("samtools index %s_pagit.bam" % name, shell=True)
    subprocess.check_call("rm 1.bam 2.bam 3.bam", shell=True)

def run_gatk(reference, processors, name, gatk):
    args = ['java', '-jar', '%s' % gatk, '-T', 'UnifiedGenotyper',
            '-R', '%s' % reference, '-nt', '%s' % processors,'-S', 'silent',
            '-mbq', '17', '-ploidy', '1', '-out_mode', 'EMIT_VARIANTS_ONLY',
            '-stand_call_conf', '100', '-stand_emit_conf', '100', '-I', '%s_pagit.bam' % name,
            '-rf', 'BadCigar']
    try:
        vcf_fh = open('gatk.out', 'w')
    except:
        print 'could not open gatk file'
    try:
        log_fh = open('gatk.log', 'w')
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
    return to_fix

def fasta_to_tab(fasta):
    infile = open(fasta, "rU")
    outfile = open("out.tab", "w")
    for record in SeqIO.parse(infile, "fasta"):
        print >> outfile, record.id, record.seq
    infile.close()
    outfile.close()

def rename_multifasta(fasta_in, prefix, fasta_out):
    """rename mutli-fasta to something meaningful"""
    handle = open(fasta_out, "w")
    for record in SeqIO.parse(open(fasta_in), "fasta"):
        print >> handle, ">"+prefix+"_"+str(autoIncrement())
        print >> handle, record.seq
    handle.close()

def fix_assembly(in_tab, to_fix, name):
    output_handle=open("%s_corrected_assembly.fasta" % name,"w")
    for line in open(in_tab, "U"):
        fields=line.split()
        seq_list=[]
        seq_list=list(fields[1])
        for x in to_fix:
            if fields[0]== x[0]:
                seq_list[int(x[1])-1]=x[2]
        print >> output_handle,">"+fields[0]
        print >> output_handle, "".join(seq_list)
    output_handle.close()

def run_loop(fileSets,error_corrector,processors,keep,coverage,proportion,start_path):
    files_and_temp_names = [(str(idx), list(f))
                            for idx, f in fileSets.iteritems()]
    lock = threading.Lock()
    def _perform_workflow(data):
        idx, f = data
        if int(get_sequence_length(f[0], idx))<=200:
	    args=['java','-jar','%s' % TRIM_PATH,'PE',
	      '%s' % f[0], '%s' % f[1], '%s.F.paired.fastq.gz' % idx, 'F.unpaired.fastq.gz',
	      '%s.R.paired.fastq.gz' % idx, 'R.unpaired.fastq.gz', '%s/bin/illumina_adapters_all.fasta:2:30:10' % UGAP_PATH,
	      'MINLEN:80']
	    try:
                vcf_fh = open('%s.trimmomatic.out' % idx, 'w')
            except:
                log_isg.logPrint('could not open trimmomatic file')
            try:
                log_fh = open('%s.trimmomatic.log' % idx, 'w')
            except:
                log_isg.logPrint('could not open log file')
	    try:
	        trim = Popen(args, stderr=vcf_fh, stdout=log_fh)
                trim.wait()
	    except:
		log_isg.logPrint("problem encountered with trimmomatic")
            """assemble sequences with spades"""
            if error_corrector=="hammer":
                subprocess.check_call("spades.py -o %s.spades -t %s -k 21,33,55,77 --careful -1 %s.F.paired.fastq.gz -2 %s.R.paired.fastq.gz  > /dev/null 2>&1" % (idx,processors,idx,idx), shell=True)
            elif error_corrector=="musket":
                ab = subprocess.call(['which', 'musket'])
                if ab == 0:
                    pass
                else:
                    print "musket isn't in your path, but needs to be!"
                    sys.exit()
                subprocess.check_call("musket -k 17 8000000 -p %s -omulti %s -inorder %s.F.paired.fastq.gz %s.R.paired.fastq.gz > /dev/null 2>&1" % (processors,idx,idx,idx), shell=True)
		subprocess.check_call("mv %s.0 %s.0.musket.fastq.gz" % (idx,idx), shell=True)
		subprocess.check_call("mv %s.1 %s.1.musket.fastq.gz" % (idx,idx), shell=True)
                subprocess.check_call("spades.py -o %s.spades -t %s -k 21,33,55,77 --only-assembler --careful -1  %s.0.musket.fastq.gz -2 %s.1.musket.fastq.gz > /dev/null 2>&1" % (idx,processors,idx,idx), shell=True)
            else:
                subprocess.check_call("spades.py -o %s.spades -t %s -k 21,33,55,77 --only-assembler --careful -1 %s.F.paired.fastq.gz -2 %s.R.paired.fastq.gz > /dev/null 2>&1" % (idx,processors,idx,idx), shell=True)
        elif int(get_sequence_length(f[0], idx))>200:
	    args=['java','-jar','%s' % TRIM_PATH,'PE',
	          '%s' % f[0], '%s' % f[1], '%s.F.paired.fastq.gz' % idx, 'F.unpaired.fastq.gz',
	          '%s.R.paired.fastq.gz' % idx, 'R.unpaired.fastq.gz', 'ILLUMINACLIP:/home/jsahl/illumina_adapters_all.fasta:2:30:10',
	          'MINLEN:150']
	    try:
                vcf_fh = open('%s.trimmomatic.out' % idx, 'w')
            except:
                log_isg.logPrint('could not open trimmomatic file')
            try:
                log_fh = open('%s.trimmomatic.log' % idx, 'w')
            except:
                log_isg.logPrint('could not open log file')
	    try:
	        trim = Popen(args, stderr=vcf_fh, stdout=log_fh)
                trim.wait()
	    except:
	        log_isg.logPrint("problem encountered with trimmomatic")
            """assemble sequences with spades"""
            if error_corrector=="hammer":
                subprocess.check_call("spades.py -o %s.spades -t %s -k 21,33,55,77,127 --careful -1 %s.F.paired.fastq.gz -2 %s.R.paired.fastq.gz  > /dev/null 2>&1" % (idx,processors,idx,idx), shell=True)
            elif error_corrector=="musket":
                ab = subprocess.call(['which', 'musket'])
                if ab == 0:
                    pass
                else:
                    print "musket isn't in your path, but needs to be!"
                    sys.exit()
                subprocess.check_call("mv %s.0 %s.0.musket.fastq.gz" % (idx,idx), shell=True)
	        subprocess.check_call("mv %s.1 %s.1.musket.fastq.gz" % (idx,idx), shell=True)
                subprocess.check_call("spades.py -o %s.spades -t %s -k 21,33,55,77,127 --only-assembler --careful -1  %s.0.musket.fastq.gz -2 %s.1.musket.fastq.gz > /dev/null 2>&1" % (idx,processors,idx,idx), shell=True)
            else:
                subprocess.check_call("spades.py -o %s.spades -t %s -k 21,33,55,77,127 --only-assembler --careful -1 %s.F.paired.fastq.gz -2 %s.R.paired.fastq.gz > /dev/null 2>&1" % (idx,processors,idx,idx), shell=True)
        else:
            pass
	os.system("zcat %s.F.paired.fastq.gz > %s_1.fastq" % (idx,idx))
        os.system("zcat %s.R.paired.fastq.gz > %s_2.fastq" % (idx,idx))
        os.system("mv %s_1.fastq %s_2.fastq %s.spades" % (idx,idx,idx))
        #os.chdir("%s.spades" % idx)
	os.system("cp %s.spades/contigs.fasta %s.spades.assembly.fasta" % (idx,idx))
        filter_seqs("%s.spades.assembly.fasta" % idx, keep, idx)
        """remove redundancies - will likely change in the near future"""
        os.system("gt sequniq -rev yes -o %s.%s.nr.spades.assembly.fasta %s.%s.spades.assembly.fasta" % (idx,keep,idx,keep))
        lock.acquire()
        run_image("%s.%s.nr.spades.assembly.fasta" % (idx,keep), "%s" % idx, get_sequence_length(f[0], idx))
        os.system("restartIMAGE.pl %s2 %s 2 partitioned > out.image2.txt 2> /dev/null" % (idx,int(get_sequence_length(f[0], idx))-30))
        os.system("restartIMAGE.pl %s4 %s 2 partitioned > out.image3.txt 2> /dev/null" % (idx,int(get_sequence_length(f[0], idx))-40))
        os.system("cp %s6/new.fa ./%s.%s.image.fasta" % (idx,idx,keep))
        lock.release()
        os.system("icorn.start.sh %s.%s.image.fasta 1 10 %s_1.fastq %s_2.fastq 200,400 300 > out.icorn.txt 2> /dev/null" % (idx,keep,idx,idx))
        #os.system("rm -rf *.gff")
	clean_fasta("%s.%s.image.fasta.11" % (idx,keep),"%s_pagit.fasta" % idx)
        subprocess.check_call("bwa index %s_pagit.fasta > /dev/null 2>&1" % idx, shell=True)
        os.system("samtools faidx %s_pagit.fasta" % idx)
        run_bwa("%s_1.fastq" % idx, "%s_2.fastq" % idx, processors, idx,"%s_pagit.fasta" % idx)
        make_bam("%s.sam" % idx, idx)
        """need to define the GATK_PATH variable"""
	os.system("java -jar %s R=%s_pagit.fasta O=%s_pagit.dict > /dev/null 2>&1" % (PICARD_PATH, idx, idx))
        run_gatk("%s_pagit.fasta" % idx, processors, idx, "%s" % GATK_PATH)
        lock.acquire()
	to_fix=parse_vcf("gatk.out", coverage, proportion)
        log_isg.logPrint("number of SNPs to fix in %s = %s" % (idx,len(to_fix)))
	if int(len(to_fix))>=1:
            fasta_to_tab("%s_pagit.fasta" % idx)
            fix_assembly("out.tab", to_fix, idx)
            try:
                rename_multifasta("%s_corrected_assembly.fasta" % idx,idx,"%s_corrected_assembly_final.fasta" % idx)
            except:
                print "error correction failed for some reason"
	    os.system("cp %s_corrected_assembly_final.fasta %s_final_assembly.fasta" % (idx,idx))
        else:
	    os.system("cp %s_pagit.fasta %s_final_assembly.fasta" % (idx,idx))
        lock.release()
	os.system("prokka --prefix %s --locustag %s --compliant --mincontiglen %s --strain %s %s_final_assembly.fasta > /dev/null 2>&1" % (idx,idx,keep,idx,idx))
        rename_multifasta("%s_final_assembly.fasta" % idx,idx,"%s_corrected_assembly_final.fasta" % idx)
	os.system("cp %s_corrected_assembly_final.fasta %s/UGAP_assembly_results" % (idx,start_path))
	os.system("cp %s/*.* %s/UGAP_assembly_results" % (idx,start_path))
	os.chdir("%s" % start_path)
    results = set(p_func.pmap(_perform_workflow,
                              files_and_temp_names,
                              num_workers=processors))
           
def main(directory,error_corrector,processors,keep,coverage,proportion,temp_files):
    """1. spades (include in bin)
       3. musket (include in bin)
       4. PAGIT (include in bin)
       5. gt (remove?)
       6. bwa (include in bin)
       7. Trimmomatic (include in bin)
       8. samtools (include in bin)
       9. GATK (include in bin)
       10. Picard tools (include in bin)"""
    start_dir = os.getcwd()
    start_path = os.path.abspath("%s" % start_dir)
    try:
        os.makedirs('%s/UGAP_assembly_results' % start_path)
        os.makedirs('%s/work_directory' % start_path)
    except OSError, e:
             if e.errno != errno.EEXIST:
                 raise 
    dir_path=os.path.abspath("%s" % directory)
    os.system("ln -s %s %s/work_directory" % (dir_path, start_path))
    fileSets=read_file_sets("%s/work_directory" % start_path)
    os.chdir("%s/work_directory" % start_path)
    log_isg.logPrint("starting loop")
    run_loop(fileSets,error_corrector,processors,keep,coverage,proportion,start_path)
    log_isg.logPrint("loop finished")
    log_isg.logPrint("cleaning up")    
    os.chdir("%s" % start_path)
    if temp_files == "F":
        os.system("rm -rf *.spades")
    else:
        pass

if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = OptionParser(usage=usage) 
    parser.add_option("-d", "--directory", dest="directory",
                      help="directory to where .fastq.gz files are found [REQUIRED]",
                      action="callback", callback=test_dir, type="string")
    parser.add_option("-e", "--error", dest="error_corrector",
                      help="error corrector, choose from musket,hammer, or none, defaults to hammer",
                      action="callback", callback=test_options, type="string", default="hammer")
    parser.add_option("-p", "--processors", dest="processors",
                      help="# of processors to use - defaults to 2",
                      default="2", type="int")
    parser.add_option("-k", "--keep", dest="keep",
                      help="minimum length of contigs to keep, defaults to 200",
                      default="200", type="int")
    parser.add_option("-c", "--coverage", dest="coverage",
                      help="minimum coverage required for correcting SNPs, defaults to 3",
                      default="3", type="int")
    parser.add_option("-i", "--proportion", dest="proportion",
                      help="minimum required proportion, defaults to 0.9",
                      action="store", type="float", default="0.9")
    parser.add_option("-t", "--temp_files", dest="temp_files",
                      help="Keep temp files? Defaults to F",
                      action="callback", callback=test_truths, type="string", default="F")
    options, args = parser.parse_args()
    
    mandatories = ["directory"]
    for m in mandatories:
        if not options.__dict__[m]:
            print "\nMust provide %s.\n" %m
            parser.print_help()
            exit(-1)

    main(options.directory,options.error_corrector,options.processors,options.keep,
         options.coverage,options.proportion,options.temp_files)

