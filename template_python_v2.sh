#!/usr/bin/env python

''' Usage:
        -i or --infi
        -o or --outfi
'''

from __future__ import print_function, division
import pysam
import getopt
import sys
import re
import collections
import operator

def getOptions(argv):
    infile      = ''
    outfile     = ''
    try:
        opts, args = getopt.getopt(argv[1:], "hi:o:", ["help", "infi=", "outfi="])
    except getopt.GetoptError:
        print( 'Usage: %s -i <infi> -o <outfi>' % argv[0], file =  sys.stderr )
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print( 'Usage: %s -i <infi> -o <outfi>' % argv[0], file = sys.stderr )
            sys.exit()
        elif opt in ("-i", "--infi"):
            infile  = arg
        elif opt in ("-o", "--outfi"):
            outfile = arg
        else:
            assert False, "unhandled option"
    return infile, outfile

# Get options:
infile, outfile = getOptions(sys.argv)

