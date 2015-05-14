#!/usr/bin/env python

"""the universal genome assembly pipeline"""

from __future__ import division
from optparse import OptionParser
import sys
import os
import re
from igs.utils import logging as log_isg
import errno
try:
    from ugap.util import *
    from igs.utils import logging as log_isg
    import ugap_single import *
except:
    print "Environment not set correctly, correct ugap_single.py environment"
    sys.exit()

"""UGAP_PATH must be modified for the install
location on your machine"""
UGAP_PATH="/Users/jsahl/UGAP"
sys.path.append('%s' % UGAP_PATH)
"""set the paths for all of the java dependencies"""
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

def parse_config_file(config_file):
    datasets = ()
    infile = open(config_file, "U")
    for line in infile:
        if line.startswith("#"):
            pass
        else:
            fields = line.split()
            datasets=((fields[0],fields[1],fields[2],fields[3],fields[4],fields[5],fields[6],fields[7],fields[8],fields[9],fields[10],fields[11],fields[12],fields[13]),)+datasets
    return datasets
          
def main(config_file, processors, memory):
    start_dir = os.getcwd()
    start_path = os.path.abspath("%s" % start_dir)
    config_path = os.path.abspath("%s" % config_file)
    try:
        os.makedirs('%s/UGAP_assembly_results' % start_path)
        os.makedirs('%s/work_directory' % start_path)
    except OSError, e:
        if e.errno != errno.EEXIST:raise 
    dir_path=os.path.abspath("%s" % directory)
    if "NULL" != reduce:
        reduce_path=os.path.abspath("%s" % reduce)
    #os.system("ln -s %s/* %s/work_directory" % (dir_path, start_path))
    #fileSets=read_file_sets("%s/work_directory" % start_path)
    os.chdir("%s/work_directory" % start_path)
    effective_jobs = int(int(memory)/16)
    if effective_jobs <=1:
        effective_jobs = 1
    effective_processors = int(int(processors)/effective_jobs)
    parameters = parse_job_parameters(config_path)
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
    parser.add_option("-c", "--config", dest="config_file",
                      help="config file that populates the UGAP single assembly",
                      action="callback", callback=test_file, type="string")
    parser.add_option("-p", "--processors", dest="processors",
                      help="# of processors available on your system, defaults to ",
                      default="2", type="int")
    parser.add_option("-m", "--memory", dest="memory",
                      help="amount of memory on the system that you using in Gb, defaults to 48(Gb)",
                      action="store", type="string", default="48")
    options, args = parser.parse_args()
    
    mandatories = ["config_file"]
    for m in mandatories:
        if not options.__dict__[m]:
            print "\nMust provide %s.\n" %m
            parser.print_help()
            exit(-1)

    main(options.config_file, options.processors, options.memory)
