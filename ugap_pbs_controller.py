#!/usr/bin/env python

"""Submits jobs to qsub,
using UGAP assembly singe"""

import sys
import os
from optparse import OptionParser
from popen2 import popen2

def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print '%s file cannot be opened' % option
        sys.exit()

def parse_config_file(config_file):
    datasets = ()
    infile = open(config_file, "U")
    for line in infile:
        fields = line.split()
        datasets=((fields[0],fields[1],fields[2],fields[3],fields[4],fields[5],fields[6],fields[7],fields[8],fields[9]),)+datasets
    return datasets

def send_jobs(datasets,my_mem):
    for data in datasets:
        output, input = popen2('qsub')
        job_name = "UGAP_%s" % data[0]
        walltime = "48:00:00"
        memory_redux = my_mem.replace("G","")
        if int(memory_redux)<48:
            my_q = "batch"
        else
            my_q = "hmem"
        processors = "nodes=1:ppn=%s" % data[9]
        command = "python /scratch/jsahl/tools/UGAP/ugap_single.py -n %s -f %s -v %s -e %s -k %s -c %s -i %s -t %s -r %s -p %s" % (data[0],data[1],data[2],data[3],data[4],data[5],data[6],data[7],data[8],data[9])
        memory = "mem=%s" % my_mem
        job_string = """#!/bin/bash
        #PBS -N %s
        #PBS -l walltime=%s
        #PBS -l %s
        #PBS -l %s
        #PBS -j oe
        #PBS -m a
        #PBS -q %s
        cd $PBS_O_WORKDIR
        module add blastall
        module add musket
        export PATH=/scratch/jsahl/tools/UGAP/bin:$PATH
        %s""" % (job_name, walltime, processors, memory, my_q, command)

        input.write(job_string)
        input.close()

        print job_string
        print output.read()
         
def main(config_file,memory):
    datasets=parse_config_file(config_file)
    my_mem = memory
    send_jobs(datasets,my_mem)
    
if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-c", "--config", dest="config_file",
                      help="config file that populates the UGAP single assembly",
                      action="callback", callback=test_file, type="string")
    parser.add_option("-m", "--memory", dest="memory",
                      help="amount of memory requested, defaults to 15G",
                      action="store", type="string", default="15000mb")
    options, args = parser.parse_args()
    mandatories = ["config_file"]
    for m in mandatories:
        if not options.__dict__[m]:
            print "\nMust provide %s.\n" %m
            parser.print_help()
            exit(-1)
    main(options.config_file,options.memory)
