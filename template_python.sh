from __future__ import print_function, division
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import GC as gc
from copy import deepcopy
from matplotlib.backends.backend_pdf import PdfPages
import glob
import re
import os
import sys
import numpy as np
import argparse
import collections
import pandas as pd
import matplotlib.pyplot as plt
# import primer3 # Since there is no primer3 for windows, this import happens after parsargs and depends on option --usePrimer3.


THIS_SCRIPT = os.path.basename(__file__)


def main(debugFlag, myInput=None, myThreshold=None):

    if debugFlag == True:
        sys.stderr.write('{0}: in debug mode\n'.format(THIS_SCRIPT))
    
    if myInput == None:
        sys.stderr.write(THIS_SCRIPT + ': --my_input_file <path/to/input> is not provided!!!\n')
        os._exit(1)
    else:
        message = '{0}: --my_input_file <{1}>'.format(THIS_SCRIPT, myInput) 
        sys.stderr.write(message + '\n')

    if myThreshold == None:
        sys.stderr.write('{0}: no threshold provided\n'.format(THIS_SCRIPT))
    else:
        message = '{0}: threshold provided <{1}>'.format(THIS_SCRIPT, myThreshold) 
        sys.stderr.write(message + '\n')

    return 0


if __name__ ==  "__main__":
 
    parser = argparse.ArgumentParser(description="Count primer and adapter dimers")
    
    # one way of setting default
    parser.add_argument('--debug', action='store_true')
    parser.set_defaults(debug=False)
    parser.add_argument('--my_input_file', help='All primer sequences in the run')
    # alternative way of setting default
    parser.add_argument('--threshold', help='Dimer output threshold', type=float, default=0.0)

    args = parser.parse_args()

    sys.exit(main(args.debug, args.my_input_file, args.threshold))
