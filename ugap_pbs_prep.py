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

def test_dir(option, opt_str, value, parser):
    if os.path.exists(value):
        setattr(parser.values, option.dest, value)
    else:
        print "directory of fastqs cannot be found"
        sys.exit()

def main(directory):
    dir_path=os.path.abspath("%s" % directory)
    fileSets=read_file_sets("%s" % dir_path)
    for k,v in fileSets.iteritems():
        print k+"\t"+'\t'.join(v)
    
if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = OptionParser(usage=usage) 
    parser.add_option("-d", "--directory", dest="directory",
                      help="directory to where .fastq.gz files are found [REQUIRED]",
                      action="callback", callback=test_dir, type="string")
    options, args = parser.parse_args()
    
    mandatories = ["directory"]
    for m in mandatories:
        if not options.__dict__[m]:
            print "\nMust provide %s.\n" %m
            parser.print_help()
            exit(-1)

    main(options.directory)
