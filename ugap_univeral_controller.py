#!/usr/bin/env python


"""one controller to rule them all"""

import sys
import os
from optparse import OptionParser
#from popen2 import popen2

def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print('%s file cannot be opened' % option)
        sys.exit()

def test_options(option, opt_str, value, parser):
    if "slurm" in value:
        setattr(parser.values, option.dest, value)
    elif "torque" in value:
        setattr(parser.values, option.dest, value)
    elif "sge" in value:
        setattr(parser.values, option.dest, value)
    else:
        print("option not supported, choose from slurm, torque, or sge")
        sys.exit()

def parse_config_file(config_file):
    datasets = ()
    with open(config_file) as infile:
        for line in infile:
            if line.startswith("#"):
                pass
            else:
                fields = line.split()
                """fields[13] is the assembler option"""
                datasets=((fields[0],fields[1],fields[2],fields[3],fields[4],fields[5],fields[6],fields[7],fields[8],fields[9],fields[10],fields[11],fields[12],fields[13],fields[14]),)+datasets
    return datasets

def send_jobs(datasets,my_mem,controller,time):
    for data in datasets:
        job_name = "UGAP_%s" % data[0]
        walltime = "%s:00:00" % time
        command = "python %s/ugap_single.py -c %s -a %s -n %s -f %s -v %s -e %s -k %s -t %s -r %s -p %s -x %s -z %s -b %s -o %s -j %s" % (data[9],data[14],data[13],data[0],data[1],data[2],data[3],data[4],data[5],data[6],data[7],data[8],data[9],data[10],data[11],data[12])
        if controller == "slurm":
            memory = "mem=%s" % my_mem
            job_string = "--job-name %s -c %s --time %s --mem %s --wrap" % (job_name,data[7],walltime,my_mem)
#SBATCH -J %s
#SBATCH -c %s
#SBATCH --time %s
#SBATCH  --mem=%s
#%s""" % (job_name,data[7],walltime,my_mem,command)

            #input.write(job_string)
            #input.close()

            print('sbatch %s "%s"' % (job_string,command))
            os.system('sbatch %s "%s"' % (job_string,command))

def main(config_file,memory,controller,time):
    datasets=parse_config_file(config_file)
    my_mem = memory
    send_jobs(datasets,my_mem, controller, time)

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
    parser.add_option("-t", "--time", dest="time",
                      help="how much time to request? Defaults to 48[h]",
                      action="store", type="int", default="48")
    options, args = parser.parse_args()
    mandatories = ["config_file"]
    for m in mandatories:
        if not options.__dict__[m]:
            print("\nMust provide %s.\n" %m)
            parser.print_help()
            exit(-1)
    main(options.config_file,options.memory,options.controller,options.time)
