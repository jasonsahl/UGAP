#!/usr/bin/env python

"""Submits jobs to sge,
using UGAP assembly single"""

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

def send_jobs(datasets):
    for data in datasets:
        output, input = popen2('qsub')
        command = "/usr/local/packages/Python-2.7/bin/python /local/projects/drasko/jsahl.dir/UGAP-master/ugap_single.py -n %s -f %s -v %s -e %s -k %s -c %s -i %s -t %s -r %s -p %s" % (data[0],data[1],data[2],data[3],data[4],data[5],data[6],data[7],data[8],data[9])
        job_name = "UGAP_%s" % data[0]
        job_string = """
#!/bin/bash
#$ -N %s
#$ -P drasko-lab
#$ -l mem_free=20G
#$ -cwd
%s""" % (job_name,command)

        input.write(job_string)
        input.close()

        print job_string
        print output.read()
 
def main(config_file):
    datasets=parse_config_file(config_file)
    send_jobs(datasets)
    
if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-c", "--config", dest="config_file",
                      help="config file that populates the UGAP single assembly",
                      action="callback", callback=test_file, type="string")
    options, args = parser.parse_args()
    mandatories = ["config_file"]
    for m in mandatories:
        if not options.__dict__[m]:
            print "\nMust provide %s.\n" %m
            parser.print_help()
            exit(-1)
    main(options.config_file)
