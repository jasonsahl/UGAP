#!/usr/bin/env python

"""ugap single, meant to be run
in conjunction with PBS, SLURM, or SGE"""

from optparse import OptionParser
import os
import sys
import subprocess
try:
    from ugap.util import *
except:
    print("Environment not set correctly, correct ugap_single.py environment")
    sys.exit()
import errno
from subprocess import Popen

def test_dir(option, opt_str, value, parser):
    if os.path.exists(value):
        setattr(parser.values, option.dest, value)
    else:
        print("directory of fastas cannot be found")
        sys.exit()

def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print('%s file cannot be opened' % option)
        sys.exit()

def test_options(option, opt_str, value, parser):
    if "hammer" in value:
        setattr(parser.values, option.dest, value)
    elif "none" in value:
        setattr(parser.values, option.dest, value)
    else:
        print("select from hammer or none")
        sys.exit()

def test_truths(option, opt_str, value, parser):
    if "T" in value:
        setattr(parser.values, option.dest, value)
    elif "F" in value:
        setattr(parser.values, option.dest, value)
    else:
        print("must select from T or F")
        sys.exit()

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

def main(forward_read,name,reverse_read,error_corrector,keep,temp_files,reduce,processors,
    careful,ugap_path,blast_nt,cov_cutoff,filter_phiX,assembler,adapter_trimmer,sample_type):
    UGAP_PATH=ugap_path
    PICARD_PATH=UGAP_PATH+"/bin/"
    TRIM_PATH=UGAP_PATH+"/bin/trimmomatic.jar"
    #updated to 1.20 on September 27th, 2016
    PILON_PATH=UGAP_PATH+"/bin/pilon-1.24.jar"
    if os.path.exists(UGAP_PATH):
        sys.path.append("%s" % UGAP_PATH)
    else:
        print("your UGAP path is not correct.  Edit the path in ugap_pbs_prep.py and try again")
        sys.exit()
    start_dir = os.getcwd()
    start_path = os.path.abspath("%s" % start_dir)
    forward_path = os.path.abspath("%s" % forward_read)
    reverse_path = os.path.abspath("%s" % reverse_read)
    blast_nt_path = os.path.abspath("%s" % blast_nt)
    if os.path.exists('%s/UGAP_assembly_results' % start_path):
        pass
    else:
        try:
            os.makedirs('%s/UGAP_assembly_results' % start_path)
        except:
            pass
    if os.path.exists('%s/%s.work_directory' % (start_path,name)):
        pass
    else:
        os.makedirs('%s/%s.work_directory' % (start_path,name))
    if "NULL" != reduce:
        reduce_path=os.path.abspath("%s" % reduce)
    """test for dependencies"""
    if "spades" in assembler:
        dependencies = ['bwa','samtools','spades.py','genomeCoverageBed','seqtk']
    else:
        dependencies = ['bwa','samtools','skesa','genomeCoverageBed','seqtk']
    if "NULL" not in blast_nt:
        rx = subprocess.call(['which', 'blastn'])
        if rx == 0:
            pass
        else:
            print("you asked to use a blastdb, but you don't have blastn in your PATH. Expect errors!")
    for dependency in dependencies:
        ra = subprocess.call(['which', '%s' % dependency])
        if ra == 0:
            pass
        else:
            print("%s is not in your path, but needs to be!" % dependency)
            sys.exit()
    """done checking for dependencies"""
    os.chdir("%s/%s.work_directory" % (start_path,name))
    #This is where I really need to focus
    if "NULL" not in reduce:
        if "NULL" not in blast_nt:
            run_single_loop(assembler,forward_path,reverse_path,name,error_corrector,processors,keep,start_path,reduce_path,careful,UGAP_PATH,TRIM_PATH,PICARD_PATH,PILON_PATH,blast_nt_path,cov_cutoff,filter_phiX,sample_type)
        else:
            run_single_loop(assembler,forward_path,reverse_path,name,error_corrector,processors,keep,start_path,reduce_path,careful,UGAP_PATH,TRIM_PATH,PICARD_PATH,PILON_PATH,blast_nt,cov_cutoff,filter_phiX,sample_type)
    else:
        if "NULL" not in blast_nt:
            run_single_loop(assembler,forward_path,reverse_path,name,error_corrector,processors,keep,start_path,reduce,careful,UGAP_PATH,TRIM_PATH,PICARD_PATH,PILON_PATH,blast_nt_path,cov_cutoff,filter_phiX,sample_type)
        else:
            run_single_loop(assembler,forward_path,reverse_path,name,error_corrector,processors,keep,start_path,reduce,careful,UGAP_PATH,TRIM_PATH,PICARD_PATH,PILON_PATH,blast_nt,cov_cutoff,filter_phiX,sample_type)
    os.chdir("%s" % start_path)
    if temp_files == "F":
        os.system("rm -rf %s.work_directory" % name)

if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-a", "--assembler", dest="assembler",
                      help="Which assembler to use? Choose from spades [default] or skesa",
                      action="store", type="string", default="spades")
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
                      action="callback", callback=test_options, type="string", default="none")
    parser.add_option("-k", "--keep", dest="keep",
                      help="minimum length of contigs to keep, defaults to 500",
                      default="500", type="int")
    parser.add_option("-t", "--temp_files", dest="temp_files",
                      help="Keep temp files? Defaults to F",
                      action="callback", callback=test_truths, type="string", default="F")
    parser.add_option("-r", "--reduce", dest="reduce",
                      help="Keep reads that don't align to provided genome",
                      action="store", type="string", default="NULL")
    parser.add_option("-p", "--processors", dest="processors",
                      help="number of processors to apply to the assembly",
                      action="store", type="int", default="4")
    parser.add_option("-x", "--careful", dest="careful",
                      help="use careful option in spades? Defaults to T",
                      action="callback", callback=test_truths, type="string", default="T")
    parser.add_option("-z", "--ugap_path", dest="ugap_path",
                      help="path to UGAP [REQUIRED]",
                      action="callback", callback=test_dir, type="string")
    parser.add_option("-b", "--blast_nt", dest="blast_nt",
                      help="PATH to blast nt database, defaults to NULL",
                      action="store", type="string", default="NULL")
    parser.add_option("-o", "--cov_cutoff", dest="cov_cutoff",
                      help="value to pass to SPAdes cov_cutoff option, default is 'auto', but can also be 'off'",
                      action="store", type="string", default="auto")
    parser.add_option("-j", "--filter_phiX", dest="filter_phiX",
                      help="use USEARCH to filter phiX from reads? defaults to T",
                      action="store", type="string", default="T")
    parser.add_option("-c", "--adapter_trimmer", dest="adapter_trimmer",
                      help="adapter trimmer, choose from trimmomatic or bbduk[default]",
                      action="store", type="string", default="T")
    parser.add_option("-d", "--sample_type", dest="sample_type",
                      help="PE or SE?",
                      action="store", type="string", default="T")
    options, args = parser.parse_args()
    mandatories = ["forward_read","name","ugap_path"]
    for m in mandatories:
        if not options.__dict__[m]:
            print("\nMust provide %s.\n" %m)
            parser.print_help()
            exit(-1)
    main(options.forward_read,options.name,options.reverse_read,options.error_corrector,options.keep,
         options.temp_files,options.reduce,options.processors,options.careful,options.ugap_path,options.blast_nt,
         options.cov_cutoff,options.filter_phiX,options.assembler,options.adapter_trimmer,options.sample_type)
