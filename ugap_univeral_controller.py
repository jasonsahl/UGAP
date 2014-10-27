#!/usr/bin/env python


"""one controller to rule them all"""

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
            datasets=((fields[0],fields[1],fields[2],fields[3],fields[4],fields[5],fields[6],fields[7],fields[8],fields[9],fields[10],fields[11],fields[12]),)+datasets
    return datasets

def send_jobs(datasets,my_mem,controller,queue):
    for data in datasets:
        if controller == "slurm":
            output, input = popen2('sbatch')
        else:
            output, input = popen2('qsub')
        job_name = "UGAP_%s" % data[0]
        walltime = "48:00:00"
        command = "python %s/ugap_single.py -n %s -f %s -v %s -e %s -k %s -c %s -i %s -t %s -r %s -p %s -x %s -z %s -b %s" % (data[11],data[0],data[1],data[2],data[3],data[4],data[5],data[6],data[7],data[8],data[9],data[10],data[11],data[12])
        if controller == "slurm":
            memory = "mem=%s" % my_mem
            job_string = \
"""#!/bin/sh
#SBATCH -p %s
#SBATCH -J %s
#SBATCH -c %s
#SBATCH  --mem=%s
%s""" % (queue, job_name, data[9], my_mem, command)

            input.write(job_string)
            input.close()

            print job_string
            print output.read()
        elif controller == "torque":
            memory = "mem=%s" % my_mem
            processors = "nodes=1:ppn=%s" % data[9]
            job_string = """#!/bin/bash
#PBS -N %s
#PBS -l walltime=%s
#PBS -l %s
#PBS -l mem=%s
#PBS -j oe
#PBS -m a
#PBS -q %s
cd $PBS_O_WORKDIR
%s""" % (job_name, walltime, processors, my_mem, queue, command)
            input.write(job_string)
            input.close()

            print job_string
            print output.read()
        elif controller == "sge":
            memory = "mem_free=%s" % my_mem
            job_string = """
#!/bin/bash
#$ -N %s
#$ -P drasko-lab
#$ -l %s
#$ -cwd
%s""" % (job_name,memory,command)

            input.write(job_string)
            input.close()

            print job_string
            print output.read()

def main(config_file, memory, controller, queue):
    datasets=parse_config_file(config_file)
    my_mem = memory
    send_jobs(datasets,my_mem, controller, queue)

if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-c", "--config", dest="config_file",
                      help="config file that populates the UGAP single assembly",
                      action="callback", callback=test_file, type="string")
    parser.add_option("-m", "--memory", dest="memory",
                      help="amount of memory requested, defaults to 15G",
                      action="store", type="string", default="15000mb")
    parser.add_option("-o", "--controller", dest="controller",
                      help="system to use: choose from slurm, torque [default], or sge",
                      action="callback", callback=test_options, type="string", default="torque")
    parser.add_option("-q", "--queue", dest="queue",
                      help="which queue to use?",
                      action="store", type="string")
    options, args = parser.parse_args()
    mandatories = ["config_file"]
    for m in mandatories:
        if not options.__dict__[m]:
            print "\nMust provide %s.\n" %m
            parser.print_help()
            exit(-1)
    main(options.config_file,options.memory,options.controller,options.queue)
        
