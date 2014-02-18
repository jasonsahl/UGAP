import os
import logging
from Bio import SeqIO
import subprocess
from subprocess import Popen
import glob
from igs.threading import functional as p_func
import threading
import decimal


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
    subprocess.check_call("samtools view -h -b -S %s -o %s.1.bam 2> /dev/null" % (in_sam, name), shell=True)
    subprocess.check_call("samtools view -u -h -F4 -o %s.2.bam %s.1.bam" % (name,name), shell=True)
    subprocess.check_call("samtools view -h -b -q1 -F4 -o %s.3.bam %s.2.bam" % (name,name), shell=True)
    subprocess.check_call("samtools sort %s.3.bam %s_renamed" % (name,name), shell=True)
    subprocess.check_call("samtools index %s_renamed.bam" % name, shell=True)
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
    return to_fix

def fasta_to_tab(fasta, idx):
    infile = open(fasta, "rU")
    outfile = open("%s.out.tab" % idx, "w")
    for record in SeqIO.parse(infile, "fasta"):
        print >> outfile, record.id, record.seq
    infile.close()
    outfile.close()

def rename_multifasta(fasta_in, prefix, fasta_out):
    """rename mutli-fasta to something meaningful"""
    rec=1
    handle = open(fasta_out, "w")
    for record in SeqIO.parse(open(fasta_in), "fasta"):
        print >> handle, ">"+prefix+"_"+str(autoIncrement())
        print >> handle, record.seq
    handle.close()

def process_coverage(name):
    curr_dir= os.getcwd()
    outfile = open("coverage_out.txt", "a")
    coverage_dict = {}
    for infile in glob.glob(os.path.join(curr_dir, "*_coverage.sample_summary")):
        for line in open(infile, "U"):
            if line.startswith("%s" % name):
                fields = line.split()
                coverage_dict.update({fields[0]:fields[2]})
    for k,v in coverage_dict.iteritems():
        print >> outfile,k,v+"\n",
                
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

def run_loop(fileSets,error_corrector,processors,keep,coverage,proportion,start_path,reduce):
    files_and_temp_names = [(str(idx), list(f))
                            for idx, f in fileSets.iteritems()]
    lock = threading.Lock()
    def _perform_workflow(data):
        idx, f = data
	if "NULL" not in reduce:
	    try:
	        subprocess.check_call("bwa index %s > /dev/null 2>&1" % reduce, shell=True)
	    except:
		print "problems with indexing input file"
		sys.exit()
	    try:
	        run_bwa("%s" % f[0], "%s" % f[1], processors, idx,"%s" % reduce)
	        os.system("samtools view -bS %s.sam > %s.bam 2> /dev/null" % (idx,idx))
	        os.system("bam2fastq -o %s#.fastq --no-aligned %s.bam > /dev/null 2>&1" % (idx,idx))
	        os.system("gzip %s_1.fastq %s_2.fastq" % (idx,idx))
	        os.system("cp %s_1.fastq.gz %s" % (idx,f[0]))
	        os.system("cp %s_2.fastq.gz %s" % (idx,f[1]))
	    except:
		print "problems depleting reads"
		sys.exit()
	else:
	    pass
        if int(get_sequence_length(f[0], idx))<=200:
            #effective_length = int(int(get_sequence_length(f[0, idx])/2))
	    args=['java','-jar','%s' % TRIM_PATH,'PE', '-threads', '%s' % processors,
	      '%s' % f[0], '%s' % f[1], '%s.F.paired.fastq.gz' % idx, 'F.unpaired.fastq.gz',
	      '%s.R.paired.fastq.gz' % idx, 'R.unpaired.fastq.gz', 'ILLUMINACLIP:%s/bin/illumina_adapters_all.fasta:2:30:10' % UGAP_PATH,
	      'MINLEN:%s' % (int(get_sequence_length(f[0],idx)/2))]
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
                try:
                    subprocess.check_call("spades.py -o %s.spades -t %s -k 21,33,55,77 --only-assembler --careful -1  %s.0.musket.fastq.gz -2 %s.1.musket.fastq.gz > /dev/null 2>&1" % (idx,processors,idx,idx), shell=True)
                except:
                    pass
            else:
                try:
                    subprocess.check_call("spades.py -o %s.spades -t %s -k 21,33,55,77 --only-assembler --careful -1 %s.F.paired.fastq.gz -2 %s.R.paired.fastq.gz > /dev/null 2>&1" % (idx,processors,idx,idx), shell=True)
                except:
                    pass
        elif int(get_sequence_length(f[0], idx))>200:
	    args=['java','-jar','%s' % TRIM_PATH,'PE',
	          '%s' % f[0], '%s' % f[1], '%s.F.paired.fastq.gz' % idx, 'F.unpaired.fastq.gz',
	          '%s.R.paired.fastq.gz' % idx, 'R.unpaired.fastq.gz', 'ILLUMINACLIP:%s/bin/illumina_adapters_all.fasta:2:30:10' % UGAP_PATH,
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
                try:
                    subprocess.check_call("spades.py -o %s.spades -t %s -k 21,33,55,77,127 --careful -1 %s.F.paired.fastq.gz -2 %s.R.paired.fastq.gz  > /dev/null 2>&1" % (idx,processors,idx,idx), shell=True)
                except:
                    pass
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
                try:
                    subprocess.check_call("spades.py -o %s.spades -t %s -k 21,33,55,77,127 --only-assembler --careful -1  %s.0.musket.fastq.gz -2 %s.1.musket.fastq.gz > /dev/null 2>&1" % (idx,processors,idx,idx), shell=True)
                except:
                    pass
            else:
                try:
                    subprocess.check_call("spades.py -o %s.spades -t %s -k 21,33,55,77,127 --only-assembler --careful -1 %s.F.paired.fastq.gz -2 %s.R.paired.fastq.gz > /dev/null 2>&1" % (idx,processors,idx,idx), shell=True)
                except:
                    pass
        else:
            pass
	try:
            os.system("zcat %s.F.paired.fastq.gz > %s_1.fastq" % (idx,idx))
            os.system("zcat %s.R.paired.fastq.gz > %s_2.fastq" % (idx,idx))
	    os.system("cp %s.spades/contigs.fasta %s.spades.assembly.fasta" % (idx,idx))
            filter_seqs("%s.spades.assembly.fasta" % idx, keep, idx)
            """remove redundancies - will likely change in the near future"""
            os.system("%s/bin/psi-cd-hit.pl -i %s.%s.spades.assembly.fasta -o %s.%s.nr.spades.assembly.fasta -c 0.99999999 -G 1 -g 1 -prog blastn -exec local -l 500" % (UGAP_PATH,idx,keep,idx,keep))
	    clean_fasta("%s.%s.nr.spades.assembly.fasta" % (idx,keep),"%s_pagit.fasta" % idx)
            rename_multifasta("%s_pagit.fasta" % idx, idx, "%s_renamed.fasta" % idx)
            subprocess.check_call("bwa index %s_renamed.fasta > /dev/null 2>&1" % idx, shell=True)
            os.system("samtools faidx %s_renamed.fasta" % idx)
            run_bwa("%s_1.fastq" % idx, "%s_2.fastq" % idx, processors, idx,"%s_renamed.fasta" % idx)
            make_bam("%s.sam" % idx, idx)
            os.system("java -jar %s/CreateSequenceDictionary.jar R=%s_renamed.fasta O=%s_renamed.dict > /dev/null 2>&1" % (PICARD_PATH, idx, idx))
            run_gatk("%s_renamed.fasta" % idx, processors, idx, "%s" % GATK_PATH)
            """run_bam_coverage stuff here"""
            os.system("java -jar %s/AddOrReplaceReadGroups.jar INPUT=%s_renamed.bam OUTPUT=%s_renamed_header.bam SORT_ORDER=coordinate RGID=%s RGLB=%s RGPL=illumina RGSM=%s RGPU=name CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT > /dev/null 2>&1" % (PICARD_PATH,idx,idx,idx,idx,idx))
            os.system("echo %s_renamed_header.bam > %s.bam.list" % (idx,idx))
            os.system("java -jar %s -R %s_renamed.fasta -T DepthOfCoverage -o %s_coverage -I %s.bam.list -rf BadCigar > /dev/null 2>&1" % (GATK_PATH,idx,idx,idx))
            process_coverage(idx)
        except:
            pass
        lock.acquire()
	try:
            to_fix=parse_vcf("%s.gatk.out" % idx, coverage, proportion)
            log_isg.logPrint("number of SNPs to fix in %s = %s" % (idx,len(to_fix)))
	    if int(len(to_fix))>=1:
                try:
                    fasta_to_tab("%s_renamed.fasta" % idx, idx)
                    fix_assembly("%s.out.tab" % idx, to_fix, idx)
                    os.system("cp %s_corrected_assembly.fasta %s_renamed.fasta" % (idx,idx))
                except:
                    print "error correction failed for some reason"
            else:
                pass

        except:
            pass
        lock.release()
        try:
            os.system("java -jar %s --genome %s_renamed.fasta --frags %s_renamed.bam --output %s_pilon > /dev/null 2>&1" % (PILON_PATH,idx,idx,idx))
	    rename_multifasta("%s_pilon.fasta" % idx, idx, "%s_final_assembly.fasta" % idx)
            os.system("prokka --prefix %s --locustag %s --compliant --mincontiglen %s --strain %s %s_final_assembly.fasta > /dev/null 2>&1" % (idx,idx,keep,idx,idx))
	    filter_seqs("%s_final_assembly.fasta" % idx, keep, idx)
            os.system("sed -i 's/\\x0//g' %s.%s.spades.assembly.fasta" % (idx,keep))
            os.system("%s/cleanFasta.pl %s.%s.spades.assembly.fasta -o %s/UGAP_assembly_results/%s_final_assembly.fasta > /dev/null 2>&1" % (PICARD_PATH,idx,keep,start_path,idx))
            os.system("cp coverage_out.txt %s/UGAP_assembly_results" % start_path)
            try:
                os.system("cp %s/*.* %s/UGAP_assembly_results" % (idx,start_path))
            except:
                pass
        except:
            pass
    results = set(p_func.pmap(_perform_workflow,
                              files_and_temp_names,
                              num_workers=processors))
