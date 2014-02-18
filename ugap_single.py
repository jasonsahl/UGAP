#!/usr/bin/env python

"""ugap single, meant to be run
in conjunction with PBS"""

from optparse import OptionParser

if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-f", "--forward", dest="forward_read",
                      help="forward read, must be *.fastq.gz",
                      action="callback", callback=test_files, type="string")
    parser.add_option("-e", "--error", dest="error_corrector",
                      help="error corrector, choose from musket,hammer, or none, defaults to hammer",
                      action="callback", callback=test_options, type="string", default="hammer")
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
    
