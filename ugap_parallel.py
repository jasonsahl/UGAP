#!/usr/bin/env python

"""the ultimate genome assembly pipeline"""

from optparse import OptionParser
import sys
import os
import re
from igs.utils import logging as log_isg
import errno
from ugap.util import *

"""UGAP_PATH must be modified for the install
location on your machine"""
UGAP_PATH="/home/jsahl/tools/UGAP"
sys.path.append('%s' % UGAP_PATH)
"""set the paths for all of the java dependencies"""
GATK_PATH=UGAP_PATH+"/bin/GenomeAnalysisTK.jar"
PICARD_PATH=UGAP_PATH+"/bin/"
TRIM_PATH=UGAP_PATH+"/bin/trimmomatic-0.30.jar"
PILON_PATH=UGAP_PATH+"/bin/pilon-1.5.jar"

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
          
def main(directory,error_corrector,processors,keep,coverage,proportion,temp_files,reduce):
    """test for dependencies"""
    start_dir = os.getcwd()
    start_path = os.path.abspath("%s" % start_dir)
    try:
        os.makedirs('%s/UGAP_assembly_results' % start_path)
        os.makedirs('%s/work_directory' % start_path)
    except OSError, e:
        if e.errno != errno.EEXIST:raise 
    try:
        os.makedirs('%s/work_directory' % start_path)
    except OSError, e:
        if e.errno != errno.EEXIST:raise 
    dir_path=os.path.abspath("%s" % directory)
    if "NULL" != reduce:
        reduce_path=os.path.abspath("%s" % reduce)
    os.system("ln -s %s/* %s/work_directory" % (dir_path, start_path))
    fileSets=read_file_sets("%s/work_directory" % start_path)
    os.chdir("%s/work_directory" % start_path)
    log_isg.logPrint("starting loop")
    if "NULL" not in reduce:
        run_loop(fileSets,error_corrector,processors,keep,coverage,proportion,start_path,reduce_path)
    else:
	run_loop(fileSets,error_corrector,processors,keep,coverage,proportion,start_path,reduce)
    log_isg.logPrint("loop finished")
    log_isg.logPrint("cleaning up")    
    os.chdir("%s" % start_path)
    if temp_files == "F":
        os.system("rm -rf *.spades")
        os.system("rm -rf work_directory")
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
    parser.add_option("-r", "--reduce", dest="reduce",
                      help="Keep reads that don't align to provided genome",
                      action="store", type="string", default="NULL")
    options, args = parser.parse_args()
    
    mandatories = ["directory"]
    for m in mandatories:
        if not options.__dict__[m]:
            print "\nMust provide %s.\n" %m
            parser.print_help()
            exit(-1)

    main(options.directory,options.error_corrector,options.processors,options.keep,
         options.coverage,options.proportion,options.temp_files,options.reduce)
