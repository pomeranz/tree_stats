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
import re

import utils.py

argparser = argparse.ArgumentParser()

argparser.add_argument("--alignment", metavar="Alignment file", type=str, required=True)
argparser.add_argument("--groups", metavar="Prefixes", type=str, required=True)
argparser.add_argument("--threads", metavar="No. of threads", type=str, required=True)
argparser.add_argument("--tree", metavar="Tree file", type=str, required=True)
#argparser.add_argument('--logdir', metavar='log_dir', type=str, required=True)
argparser.add_argument('--outroot', metavar='out_root', type=str, required=True)
#argparser.add_argument('--rerun', action='store_true')



args = argparser.parse_args()

tdg09_cmd = "java -cp dist/tdg09.jar tdg09.Analyse -alignment {0} -tree {1} -groups {3} -threads {4} > {5}"

alignment = args.alignment
groups = args.groups
threads = args.threads
tree = args.tree
#logroot = args.logdir
outroot = args.outroot

utils.check_dir(logroot)
utils.check_dir(outroot)


for infile in glob(path.join(inroot, "*", "*_prank.best.fas")):
        print infile
        basename = path.basename(infile).partition('.')[0]
        basename = "".join(basename.split("_")[0] + "_" + basename.split("_")[1])
        prefix = basename.partition('_')[0][:2]
        
        outdir = path.join(outroot, prefix)
        utils.check_dir(outdir)
        outfile = path.join(outdir, basename + ".txt")
        
        treedir = path.join(treeroot, prefix)
        treefile = path.join(treedir, basename + '.nh')
        
        tdg09 = tdg09_cmd.format(infile, treefile, groups, threads, outroot)        
        
        p = Popen(['bsub', '-R', 'rusage[tmp=512]',
               '-cwd /tmp', tdg09])
        
        p.wait()
        
        