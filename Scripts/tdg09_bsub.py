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
import glob
import os
from os import path
import argparse
import re

import utils

argparser = argparse.ArgumentParser()

argparser.add_argument("--alignment", -metavar="Alignment file", type=str, required=True)
argparser.add_argument("--groups", metavar="Prefixes", type=str, required=True)
argparser.add_argument("--threads", metavar="No. of threads", type=str, required=True)
argparser.add_argument("--tree", metavar="Tree file", type=str, required=True)
argparser.add_argument('--logdir', metavar='log_dir', type=str, required=True)
argparser.add_argument('--outroot', metavar='out_root', type=str, required=True)
argparser.add_argument('--rerun', action='store_true')



args = argparser.parse_args()

tdg09_cmd = "java -cp dist/tdg09.jar tdg09.Analyse -alignment {0} -groups {1} -threads {3} -tree {4}"

alignment = args.alignment
groups = args.groups
threads = args.threads
tree = args.tree
logroot = args.logdir
outroot = args.outroot

utils.check_dir(logroot)
utils.check_dir(outroot)

"""
GET TDG09 TO WORK ON YOUR PC AND THEN FINISH THIS SCRIPT """
