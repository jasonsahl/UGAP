#!/usr/bin/env python

"""ugap single, meant to be run
in conjunction with PBS"""

from optparse import OptionParser
import os
import sys
import subprocess
try:
    from ugap.util import *
    from igs.utils import logging as log_isg
except:
    print "Environment not set correctly"
    sys.exit()
import errno
from subprocess import Popen


UGAP_PATH="/scratch/jsahl/tools/UGAP"
sys.path.append('%s' % UGAP_PATH)
sys.path.append('%s/share' % UGAP_PATH)
GATK_PATH=UGAP_PATH+"/bin/GenomeAnalysisTK.jar"
PICARD_PATH=UGAP_PATH+"/bin/"
TRIM_PATH=UGAP_PATH+"/bin/trimmomatic-0.30.jar"
PILON_PATH=UGAP_PATH+"/bin/pilon-1.5.jar"

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

def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print '%s file cannot be opened' % option
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

def run_single_loop(forward_path,reverse_path,name,error_corrector,processors,keep,coverage,proportion,start_path,reduce):
    if "NULL" not in reduce:
        try:
            subprocess.check_call("bwa index %s > /dev/null 2>&1" % reduce, shell=True)
        except:
            print "problems with indexing input file"
            sys.exit()
        try:
            run_bwa("%s" % forward_path, "%s" % reverse_path, processors, name,"%s" % reduce)
            os.system("samtools view -bS %s.sam > %s.bam 2> /dev/null" % (name,name))
            os.system("bam2fastq -o %s#.fastq --no-aligned %s.bam > /dev/null 2>&1" % (name,name))
            os.system("gzip %s_1.fastq %s_2.fastq" % (name,name))
            os.system("cp %s_1.fastq.gz %s" % (name,forward_path))
            os.system("cp %s_2.fastq.gz %s" % (name,reverse_path))
        except:
            print "problems depleting reads"
            sys.exit()
    else:
        pass
    if int(get_sequence_length(forward_path, name))<=200:
        args=['java','-jar','%s' % TRIM_PATH,'PE', '-threads', '%s' % processors,
              '%s' % forward_path, '%s' % reverse_path, '%s.F.paired.fastq.gz' % name, 'F.unpaired.fastq.gz',
	      '%s.R.paired.fastq.gz' % name, 'R.unpaired.fastq.gz', 'ILLUMINACLIP:%s/bin/illumina_adapters_all.fasta:2:30:10' % UGAP_PATH,
	      'MINLEN:%s' % (int(get_sequence_length(forward_path,name)/2))]
        try:
            vcf_fh = open('%s.trimmomatic.out' % name, 'w')
        except:
            log_isg.logPrint('could not open trimmomatic file')
        try:
            log_fh = open('%s.trimmomatic.log' % name, 'w')
        except:
            log_isg.logPrint('could not open log file')
        try:
            trim = Popen(args, stderr=vcf_fh, stdout=log_fh)
            trim.wait()
        except:
            log_isg.logPrint("problem encountered with trimmomatic")
        if error_corrector=="hammer":
            subprocess.check_call("spades.py -o %s.spades -t %s -k 21,33,55,77 --careful -1 %s.F.paired.fastq.gz -2 %s.R.paired.fastq.gz  > /dev/null 2>&1" % (name,processors,name,name), shell=True)
        elif error_corrector=="musket":
            ab = subprocess.call(['which', 'musket'])
            if ab == 0:
                pass
            else:
                print "musket isn't in your path, but needs to be!"
                sys.exit()
            subprocess.check_call("musket -k 17 8000000 -p %s -omulti %s -inorder %s.F.paired.fastq.gz %s.R.paired.fastq.gz > /dev/null 2>&1" % (processors,name,name,name), shell=True)
            subprocess.check_call("mv %s.0 %s.0.musket.fastq.gz" % (name,name), shell=True)
            subprocess.check_call("mv %s.1 %s.1.musket.fastq.gz" % (name,name), shell=True)
            subprocess.check_call("spades.py -o %s.spades -t %s -k 21,33,55,77 --only-assembler --careful -1  %s.0.musket.fastq.gz -2 %s.1.musket.fastq.gz > /dev/null 2>&1" % (name,processors,name,name), shell=True)
        else:
            subprocess.check_call("spades.py -o %s.spades -t %s -k 21,33,55,77 --only-assembler --careful -1 %s.F.paired.fastq.gz -2 %s.R.paired.fastq.gz > /dev/null 2>&1" % (name,processors,name,name), shell=True)
    if int(get_sequence_length(forward_path, name))>200:
        args=['java','-jar','%s' % TRIM_PATH,'PE',
              '%s' % forward_path, '%s' % reverse_path, '%s.F.paired.fastq.gz' % name, 'F.unpaired.fastq.gz',
	      '%s.R.paired.fastq.gz' % name, 'R.unpaired.fastq.gz', 'ILLUMINACLIP:%s/bin/illumina_adapters_all.fasta:2:30:10' % UGAP_PATH,
	      'MINLEN:150']
        try:
            vcf_fh = open('%s.trimmomatic.out' % name, 'w')
        except:
            log_isg.logPrint('could not open trimmomatic file')
        try:
            log_fh = open('%s.trimmomatic.log' % name, 'w')
        except:
            log_isg.logPrint('could not open log file')
        try:
            trim = Popen(args, stderr=vcf_fh, stdout=log_fh)
            trim.wait()
        except:
            log_isg.logPrint("problem encountered with trimmomatic")
        """assemble sequences with spades"""
        if error_corrector=="hammer":
            subprocess.check_call("spades.py -o %s.spades -t %s -k 21,33,55,77,127 --careful -1 %s.F.paired.fastq.gz -2 %s.R.paired.fastq.gz  > /dev/null 2>&1" % (name,processors,name,name), shell=True)
        elif error_corrector=="musket":
            ab = subprocess.call(['which', 'musket'])
            if ab == 0:
                pass
            else:
                print "musket isn't in your path, but needs to be!"
                sys.exit()
            subprocess.check_call("musket -k 17 8000000 -p %s -omulti %s -inorder %s.F.paired.fastq.gz %s.R.paired.fastq.gz > /dev/null 2>&1" % (processors,name,name,name), shell=True)
            subprocess.check_call("mv %s.0 %s.0.musket.fastq.gz" % (name,name), shell=True)
            subprocess.check_call("mv %s.1 %s.1.musket.fastq.gz" % (name,name), shell=True)
            subprocess.check_call("spades.py -o %s.spades -t %s -k 21,33,55,77,127 --only-assembler --careful -1  %s.0.musket.fastq.gz -2 %s.1.musket.fastq.gz > /dev/null 2>&1" % (name,processors,name,name), shell=True)
        else:
            subprocess.check_call("spades.py -o %s.spades -t %s -k 21,33,55,77,127 --only-assembler --careful -1 %s.F.paired.fastq.gz -2 %s.R.paired.fastq.gz > /dev/null 2>&1" % (name,processors,name,name), shell=True)
    os.system("gzip -dc %s.F.paired.fastq.gz > %s_1.fastq" % (name,name))
    os.system("gzip -dc %s.R.paired.fastq.gz > %s_2.fastq" % (name,name))
    os.system("cp %s.spades/contigs.fasta %s.spades.assembly.fasta" % (name,name))
    filter_seqs("%s.spades.assembly.fasta" % name, keep, name)
    os.system("%s/bin/psi-cd-hit.pl -i %s.%s.spades.assembly.fasta -o %s.%s.nr.spades.assembly.fasta -c 0.99999999 -G 1 -g 1 -prog blastn -exec local -l 500" % (UGAP_PATH,name,keep,name,keep))
    clean_fasta("%s.%s.nr.spades.assembly.fasta" % (name,keep),"%s_pagit.fasta" % name)
    rename_multifasta("%s_pagit.fasta" % name, name, "%s_renamed.fasta" % name)
    subprocess.check_call("bwa index %s_renamed.fasta > /dev/null 2>&1" % name, shell=True)
    os.system("samtools faidx %s_renamed.fasta" % name)
    run_bwa("%s_1.fastq" % name, "%s_2.fastq" % name, processors, name,"%s_renamed.fasta" % name)
    make_bam("%s.sam" % name, name)
    os.system("java -jar %s/CreateSequenceDictionary.jar R=%s_renamed.fasta O=%s_renamed.dict > /dev/null 2>&1" % (PICARD_PATH, name, name))
    run_gatk("%s_renamed.fasta" % name, processors, name, "%s" % GATK_PATH)
    """run_bam_coverage stuff here"""
    os.system("java -jar %s/AddOrReplaceReadGroups.jar INPUT=%s_renamed.bam OUTPUT=%s_renamed_header.bam SORT_ORDER=coordinate RGID=%s RGLB=%s RGPL=illumina RGSM=%s RGPU=name CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT > /dev/null 2>&1" % (PICARD_PATH,name,name,name,name,name))
    os.system("echo %s_renamed_header.bam > %s.bam.list" % (name,name))
    os.system("java -jar %s -R %s_renamed.fasta -T DepthOfCoverage -o %s_coverage -I %s.bam.list -rf BadCigar > /dev/null 2>&1" % (GATK_PATH,name,name,name))
    process_coverage(name)
    try:
        to_fix=parse_vcf("%s.gatk.out" % name, coverage, proportion)
        log_isg.logPrint("number of SNPs to fix in %s = %s" % (name,len(to_fix)))
	if int(len(to_fix))>=1:
            try:
                fasta_to_tab("%s_renamed.fasta" % name, name)
                fix_assembly("%s.out.tab" % name, to_fix, name)
                os.system("cp %s_corrected_assembly.fasta %s_renamed.fasta" % (name,name))
            except:
                print "error correction failed for some reason"
        else:
            pass
    except:
        pass
    try:
        os.system("java -jar %s --genome %s_renamed.fasta --frags %s_renamed.bam --output %s_pilon > /dev/null 2>&1" % (PILON_PATH,name,name,name))
	rename_multifasta("%s_pilon.fasta" % name, name, "%s_final_assembly.fasta" % name)
        os.system("prokka --prefix %s --locustag %s --compliant --mincontiglen %s --strain %s %s_final_assembly.fasta > /dev/null 2>&1" % (name,name,keep,name,name))
	filter_seqs("%s_final_assembly.fasta" % name, keep, name)
        os.system("sed -i 's/\\x0//g' %s.%s.spades.assembly.fasta" % (name,keep))
        os.system("%s/cleanFasta.pl %s.%s.spades.assembly.fasta -o %s/UGAP_assembly_results/%s_final_assembly.fasta > /dev/null 2>&1" % (PICARD_PATH,name,keep,start_path,name))
        os.system("cp coverage_out.txt %s/UGAP_assembly_results/%s_coverage.txt" % (start_path,name))
        try:
            os.system("cp %s/*.* %s/UGAP_assembly_results" % (name,start_path))
        except:
            pass
    except:
        pass

def main(forward_read,name,reverse_read,error_corrector,keep,coverage,proportion,temp_files,reduce,processors):
    start_dir = os.getcwd()
    start_path = os.path.abspath("%s" % start_dir)
    forward_path = os.path.abspath("%s" % forward_read)
    reverse_path = os.path.abspath("%s" % reverse_read)
    try:
        os.makedirs('%s/UGAP_assembly_results' % start_path)
    except OSError, e:
        if e.errno != errno.EEXIST:raise
    try:
        os.makedirs('%s/%s.work_directory' % (start_path,name))
    except OSError, e:
        if e.errno != errno.EEXIST:raise
    if "NULL" != reduce:
        reduce_path=os.path.abspath("%s" % reduce)
    os.chdir("%s/%s.work_directory" % (start_path,name))
    if "NULL" not in reduce:
        run_single_loop(forward_path,reverse_path,name,error_corrector,processors,keep,coverage,proportion,start_path,reduce_path)
    else:
	run_single_loop(forward_path,reverse_path,name,error_corrector,processors,keep,coverage,proportion,start_path,reduce)
    os.chdir("%s" % start_path)
    if temp_files == "F":
        os.system("rm -rf %s.work_directory" % name)

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
                      help="reverse read, must be *.fastq.gz [REQUIRED]",
                      action="callback", callback=test_file, type="string")
    parser.add_option("-e", "--error", dest="error_corrector",
                      help="error corrector, choose from musket,hammer, or none, defaults to hammer",
                      action="callback", callback=test_options, type="string", default="hammer")
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
    parser.add_option("-r", "--reduce", dest="reduce",
                      help="Keep reads that don't align to provided genome",
                      action="store", type="string", default="NULL")
    parser.add_option("-p", "--processors", dest="processors",
                      help="number of processors to apply to the assembly",
                      action="store", type="int", default="4")
    options, args = parser.parse_args()
    mandatories = ["forward_read","name","reverse_read"]
    for m in mandatories:
        if not options.__dict__[m]:
            print "\nMust provide %s.\n" %m
            parser.print_help()
            exit(-1)
    main(options.forward_read,options.name,options.reverse_read,options.error_corrector,options.keep,options.coverage,
         options.proportion,options.temp_files,options.reduce,options.processors)
    
