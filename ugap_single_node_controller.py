#!/usr/bin/env python


"""one controller to rule them all"""

import sys
import os
from optparse import OptionParser
from popen2 import popen2
import errno

from ugap.util import run_single_loop

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
        #os.makedirs('%s/ugap_work_directory' % start_path)
    except OSError, e:
        if e.errno != errno.EEXIST:raise
    try:
        os.makedirs('%s/ugap_work_directory' % start_path)
    except OSError, e:
        if e.errno != errno.EEXIST:raise
    if "NULL" != reduce:
        reduce_path=os.path.abspath("%s" % reduce)
    effective_jobs = int(int(memory)/8000)
    if effective_jobs <=1:
        effective_jobs = 1
    effective_processors = int(int(processors)/effective_jobs)
    os.chdir("%s/ugap_work_directory" % start_dir) 
    keep_stuff = []
    def _perform_workflow(data):
        f = data
        run_single_loop(f[1],f[2],f[0],f[3],f[7],f[4],start_path,f[6],f[8],UGAP_PATH,TRIM_PATH,PICARD_PATH,PILON_PATH,f[10],f[11])
        keep_stuff.append(f[5])
    results = set(p_func.pmap(_perform_workflow,
                              datasets,
                              num_workers=effective_jobs))
    if "F" in keep_stuff:
        pass
    else:
        os.sytem("rm -rf ugap_work_directory")

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
        
