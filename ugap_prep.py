#!/usr/bin/env python

"""UGAP with PBS"""

"""Step1 - generate input files
that can then be distributed to the PBS system"""

from optparse import OptionParser
import re
import glob
import os
import logging
from ugap.util import read_file_sets
from ugap.util import get_readFile_components
import sys
import subprocess

UGAP_PATH="/Users/jasonsahl/tools/UGAP"

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

def main(directory,error_corrector,keep,temp_files,reduce,processors,careful,blast_nt,cov_cutoff,phiX_filter):
    dir_path=os.path.abspath("%s" % directory)
    fileSets=read_file_sets("%s" % dir_path)
    reduce_path = os.path.abspath("%s" % reduce)
    if phiX_filter == "T":
        dependencies = ['bwa','samtools','spades.py','usearch']
    else:
        dependencies = ['bwa','samtools','spades.py']
    for dependency in dependencies:
        ra = subprocess.check_call('which %s > /dev/null 2>&1' % dependency, shell=True)
        if ra == 0:
            pass
        else:
            print "%s is not in your path, but needs to be!" % dependency
            sys.exit()
    for k,v in fileSets.iteritems():
        print k+"\t"+'\t'.join(v)+"\t"+str(error_corrector)+"\t"+str(keep)+"\t"+str(temp_files)+"\t"+str(reduce_path)+"\t"+str(processors)+"\t"+str(careful)+"\t"+str(UGAP_PATH)+"\t"+str(blast_nt)+"\t"+str(cov_cutoff)+"\t"+str(phiX_filter)

if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-d", "--directory", dest="directory",
                      help="directory to where .fastq.gz files are found [REQUIRED]",
                      action="callback", callback=test_dir, type="string")
    parser.add_option("-e", "--error", dest="error_corrector",
                      help="error corrector, choose from musket,hammer, or none, defaults to hammer",
                      action="callback", callback=test_options, type="string", default="hammer")
    parser.add_option("-k", "--keep", dest="keep",
                      help="minimum length of contigs to keep, defaults to 200",
                      default="200", type="int")
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
                      help="PATH to BLAST nt database, defaults to NULL",
                      action="store", type="string", default="NULL")
    parser.add_option("-o", "--cov_cutoff", dest="cov_cutoff",
                      help="cov_cutoff value in SPAdes, can be integer or 'off', defaults to 'auto'",
                      action="store", type="string", default="auto")
    parser.add_option("-z", "--phiX", dest="phiX_filter",
                      help="filter for PhiX? Defaults to T, choose from T or F",
                      action="callback", callback=test_truths, type="string", default="T")
    options, args = parser.parse_args()

    mandatories = ["directory"]
    for m in mandatories:
        if not options.__dict__[m]:
            print "\nMust provide %s.\n" %m
            parser.print_help()
            exit(-1)

    main(options.directory,options.error_corrector,options.keep,options.temp_files,
         options.reduce,options.processors,options.careful,options.blast_nt,options.cov_cutoff,
         options.phiX_filter)
