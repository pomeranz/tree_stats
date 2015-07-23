# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 14:22:31 2015

@author: gideon

Description:
this script submits tdg09 jobs to the EBI cluster.

Usage:
python Scripts/tdg09_bsub.py --alignment [] --groups [] --threads [] --tree [] --logdir []

"""

from subprocess import Popen
from glob import glob
import os
from os import path
import argparse
#import re

from utils import check_dir

argparser = argparse.ArgumentParser()

argparser.add_argument("--alignment", metavar="Alignment file", type=str, required=True)
argparser.add_argument("--groups", metavar="Prefixes", type=str, required=True)
argparser.add_argument("--threads", metavar="No. of threads", type=str, required=True)
argparser.add_argument("--tree", metavar="Tree file", type=str, required=True)
argparser.add_argument('--logdir', metavar='log_dir', type=str, required=True)
argparser.add_argument('--outroot', metavar='out_root', type=str, required=True)
#argparser.add_argument('--rerun', action='store_true')



args = argparser.parse_args()

tdg09_cmd = "'java -cp Software/tdg09-1.1.1/dist/tdg09.jar tdg09.Analyse -alignment {0} -tree {1} -groups {2} -threads {3} > {4}'"

alignment = args.alignment
groups = args.groups
threads = args.threads
treeroot = args.tree
logroot = args.logdir
outroot = args.outroot

check_dir(outroot)

for infile in glob(path.join(alignment, "*", "*_prank.best.fas")):
        print infile
        basename = path.basename(infile).partition('.')[0]
        basename = "".join(basename.split("_")[0] + "_" + basename.split("_")[1])
        prefix = basename.partition('_')[0][:2]
        
        outdir = path.join(outroot, prefix)
        check_dir(outdir)
        outfile = path.join(outdir, basename + ".txt")
        
        logdir = path.join(logroot, prefix)
        check_dir(logroot)
        
        logfile = path.join(logdir, basename + '.log')
        
        
        treedir = path.join(treeroot, prefix)
        treefile = path.join(treedir, basename + '.nh')
        
        tdg09 = tdg09_cmd.format(infile, treefile, groups, threads, outroot)        
        
        p = Popen(['bsub', '-R', 'rusage[tmp=512]', '-o'+logfile,
               '-cwd /tmp', tdg09])
        
        p.wait()
        
        