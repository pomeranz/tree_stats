from glob import glob
import os
from os import path
import itertools
import re
from Bio import AlignIO 
import sys
import copy
import argparse

argparser = argparse.ArgumentParser()

argparser.add_argument('--clade', metavar='clade', type=str, required=True)
argparser.add_argument('--inroot', metavar='in_root', type=str, required=True)
argparser.add_argument('--outfile', metavar='out_file', type=str, required=True)

args = argparser.parse_args()

outfile = open(args.outfile, 'w')

colnames =  [ "ID", "pos", "alpha", "beta", "diff", "Palphasmaller", "Palphagreater", "BayesFactor", "PSRF", "Neff", "Nseqs" ]
print >>outfile, "\t".join(colnames)
for f in glob(path.join(args.inroot, args.clade, '*', '*.fubar.csv')):
    print "Processing", f
    dirname, basename = path.split(f)
    name_core = basename.partition('.')[0]
    aln = AlignIO.read(path.join(dirname, name_core + '.fa'), 'fasta')
    cols = []
    nongap = 0
    for i in range(0, aln.get_alignment_length(), 3):
        col = aln[i:(i+3), i]
        for it in col:
            if it != '---':
               nongap += 1 
                
        cols.append(str(nongap))

    infile = open(f)
    infile.readline() # Header
    for l in infile:
        f = l.rstrip().split(',')
        pos = int(f[0])
        print >>outfile, name_core + '\t' + '\t'.join(f) + '\t' + cols[pos-1]
