#!/usr/bin/env python

"""Prep UGAP for running with the controller"""

"""Step1 - generate input files
that can then be distributed to the job management system"""

from optparse import OptionParser
import re
import glob
import os
import logging
from ugap.util import read_file_sets
from ugap.util import get_readFile_components
import sys
import subprocess

UGAP_PATH="/scratch/js2829/UGAP"

def test_dir(option, opt_str, value, parser):
    if os.path.exists(value):
        setattr(parser.values, option.dest, value)
    else:
        print("directory of fastqs cannot be found")
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

def test_assembler(option, opt_str, value, parser):
    if "spades" in value:
        setattr(parser.values, option.dest, value)
    elif "skesa" in value:
        setattr(parser.values, option.dest, value)
    else:
        print("must select from skesa or spades")
        sys.exit()

def is_tool(name):
    from distutils.spawn import find_executable
    return find_executable(name) is not None

def main(assembler,directory,error_corrector,keep,temp_files,reduce,processors,careful,blast_nt,cov_cutoff,phiX_filter,adapter_trimmer):
    dir_path=os.path.abspath("%s" % directory)
    fileSets=read_file_sets("%s" % dir_path)
    reduce_path = os.path.abspath("%s" % reduce)
    #TODO:Need to change the default aligner to minimap
    dependencies = ['bwa','samtools','seqtk','polypolish']
    if phiX_filter == "T":
        dependencies.append("usearch")
    if assembler == "spades":
        dependencies.append("spades.py")
    elif assembler == "skesa":
        dependencies.append("skesa")
    if adapter_trimmer == "trimmomatic":
        dependencies.append("trimmomatic")
    elif adapter_trimmer == "bbduk":
        dependencies.append("bbduk.sh")
    for dependency in dependencies:
        result = is_tool(dependency)
        if result is False:
            print("%s isn't in your path, but needs to be!" % result)
            sys.exit()
        else:
            pass
    for k,v in fileSets.items():
        if len(v) == 2:
            print(k+"\t"+'\t'.join(v)+"\t"+str(error_corrector)+"\t"+str(keep)+"\t"+str(temp_files)+"\t"+str(reduce_path)+"\t"+str(processors)+"\t"+str(careful)+"\t"+str(UGAP_PATH)+"\t"+str(blast_nt)+"\t"+str(cov_cutoff)+"\t"+str(phiX_filter)+"\t"+str(assembler)+"\t"+str(adapter_trimmer)+"\t"+"PE")
        elif len(v) == 1:
            print(k+"\t"+'\t'.join(v)+"\t"+str(error_corrector)+"\t"+str(keep)+"\t"+str(temp_files)+"\t"+str(reduce_path)+"\t"+str(processors)+"\t"+str(careful)+"\t"+str(UGAP_PATH)+"\t"+str(blast_nt)+"\t"+str(cov_cutoff)+"\t"+str(phiX_filter)+"\t"+str(assembler)+"\t"+str(adapter_trimmer)+"\t"+"SE")

if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-a", "--assembler", dest="assembler",
                      help="which assembler to use, choose from skesa or spades [default]",
                      action="callback", callback=test_assembler, type="string", default="spades")
    parser.add_option("-d", "--directory", dest="directory",
                      help="directory to where .fastq.gz files are found [REQUIRED]",
                      action="callback", callback=test_dir, type="string")
    parser.add_option("-e", "--error", dest="error_corrector",
                      help="error corrector, choose from hammer or none, defaults to hammer",
                      action="callback", callback=test_options, type="string", default="none")
    parser.add_option("-k", "--keep", dest="keep",
                      help="minimum length of contigs to keep, defaults to 500",
                      default="500", type="int")
    parser.add_option("-t", "--temp_files", dest="temp_files",
                      help="Keep temp files? Defaults to F",
                      action="callback", callback=test_truths, type="string", default="F")
    parser.add_option("-r", "--reduce", dest="reduce",
                      help="Keep reads that don't align to provided genome, defaults to NULL",
                      action="store", type="string", default="NULL")
    parser.add_option("-p", "--processors", dest="processors",
                      help="number of processors to apply to the assembly, defaults to 4",
                      action="store", type="int", default="4")
    parser.add_option("-x", "--careful", dest="careful",
                      help="use careful option in spades? Defaults to T",
                      action="callback", callback=test_truths, type="string", default="T")
    parser.add_option("-b", "--blast_nt", dest="blast_nt",
                      help="PATH to BLAST nt database, defaults to NULL, if selected then NASP will be run",
                      action="store", type="string", default="NULL")
    parser.add_option("-o", "--cov_cutoff", dest="cov_cutoff",
                      help="cov_cutoff value in SPAdes, can be integer or 'off', defaults to 'auto'",
                      action="store", type="string", default="auto")
    parser.add_option("-z", "--phiX", dest="phiX_filter",
                      help="filter for PhiX? Defaults to T, choose from T or F",
                      action="callback", callback=test_truths, type="string", default="T")
    parser.add_option("-c", "--adapter_trimmer", dest="adapter_trimmer",
                      help="which trimmer to use: trimmomatic or bbduk [default]",
                      action="store", type="string", default="bbduk")
    options, args = parser.parse_args()

    mandatories = ["directory"]
    for m in mandatories:
        if not options.__dict__[m]:
            print("\nMust provide %s.\n" %m)
            parser.print_help()
            exit(-1)

    main(options.assembler,options.directory,options.error_corrector,options.keep,options.temp_files,
         options.reduce,options.processors,options.careful,options.blast_nt,options.cov_cutoff,
         options.phiX_filter,options.adapter_trimmer)
