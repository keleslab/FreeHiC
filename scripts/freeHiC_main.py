import pstats, cProfile
import pyximport;pyximport.install()
from freeHiC_cython import *
import argparse
import os
import sys
import re

def get_args():
    '''get arguments'''
    parser = argparse.ArgumentParser(description = '------------Usage Start------------',
                                     epilog = '------------Usage End------------')
    parser.add_argument('-f', '--fragment', help = 'Restriction Enzyme(RE) fragment file.', default = None)
    parser.add_argument('-i', '--interaction', help = 'Fragment interactions file.', default = None)
    parser.add_argument('-o', '--outdir', help = 'Output directory to save category results.', default = None)
    parser.add_argument('-fn', '--fileName', help = 'Output file name. Default will be "freeHiC.simu". ', default = "freeHiC.simu")
    parser.add_argument('-n', '--number', help = 'Total number of interactions to simulate. Default value is the number of fragment interactions in the file provided by user through -i.', default = None)
    parser.add_argument('-m', '--mutation', help = 'Percentage of mutation (0~100). Default value will be the percentage of mutation in the provided sample data.', default = 404)
    parser.add_argument('-idl', '--indel', help = 'Percentage of insertion and deletion (0~100). Default value will be the percentage of indel in the provided sample data', default = 404)
    parser.add_argument('-c', '--chimeric', help = 'Percentage of chimeric reads (0~100). Default value will be the percentage of chimeric in the provided sample data', default = 404)
    parser.add_argument('-r', '--readLen', help = 'The length of reads to be simulated.', default = None)
    parser.add_argument('-d', '--distance', help = 'Maximum distance allowed to simulate interaction coordinate away from the cloest restriction enzyme cutting site. Default is 500', default = 500)
    parser.add_argument('-s', '--bedtools', help = 'Path to bedtools executable file.', default = None)
    parser.add_argument('-g', '--genome', help = 'Path to reference genome fasta file.', default = None)
    parser.add_argument('-sf', '--summaryFile', help = '(Optional) Summary file path . Default is infileName.alignSummary.', default = None)
    parser.add_argument('-v', '--verbose', help = '(Optional) Verbose. Default is true.', default = True)

    args = parser.parse_args()

    
    if args.fragment is None or args.interaction is None or args.bedtools is None or args.genome is None or args.readLen is None:
        parser.print_help()
        sys.exit()
    if not os.path.exists(args.fragment):
        print("Resctriction Enzyme fragment file does not exist!")
        sys.exit()
    if not os.path.exists(args.interaction):
        print("Fragment interaction file does not exist!")
        sys.exit()
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    if int(args.mutation) != 404 and (int(args.mutation) > 100 or int(args.mutation) < 0):
        print("Percentage of mutation should be within [0, 100]!")
        sys.exit()
    if int(args.mutation) != 404 and (int(args.indel) > 100 or int(args.indel) < 0):
        print("Percentage of indel should be within [0, 100]!")
        sys.exit()
    if int(args.chimeric) != 404 and (int(args.chimeric) > 100 or int(args.chimeric) < 0):
        print("Percentage of chimeric reads should be within [0, 100]!")
        sys.exit()
    
    return args

if __name__ == "__main__":
    args = get_args()
    if args.verbose:
        print("RE fragment = ", args.fragment)
        print("Fragment interaction file = ", args.interaction)
        print("Output directory = ", args.outdir)
        print("Output file name = ", args.fileName)
        print("Total number of interactions to simulate = ", args.number)
        print("Read length to simulate = ", args.readLen)
        print("Maximum distance away from restriction site = ", args.distance)
        print("Path to bedtools = ", args.bedtools)
        print("Path to reference genome fasta file = ", args.genome)
        print("Path to summary file = ", args.summaryFile)
        print("Verbose = ", args.verbose)

    fragment = args.fragment
    interaction = args.interaction
    outdir = args.outdir
    fileName = args.fileName
    if args.number is None:
        number = 0
    else:
        number = int(args.number)
    mutationP = int(args.mutation)
    indelP = int(args.indel)
    chimericP = int(args.chimeric)
    readLen = int(args.readLen)
    distance = int(args.distance)
    bedtools = args.bedtools
    genome = args.genome
    summaryFile = args.summaryFile
    verbose = args.verbose
    main(fragment, interaction, outdir, fileName, number, mutationP, indelP, chimericP, readLen, distance, bedtools, genome, summaryFile, verbose)
