# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 16:37:55 2015

@author: gideon

Description:
run within prepare_slr output folder
"""
import os
from os import path
from glob import glob
import argparse
from subprocess import Popen

from utils import check_dir

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# CMD interface with argparse

argparser = argparse.ArgumentParser()

argparser.add_argument("--ctlroot", metavar="Ctl file Dir", type=str, required=True)
argparser.add_argument("--logroot", metavar="Output Dir + working dir + logdir", type=str, required=True)

args = argparser.parse_args()

ctlroot = args.ctlroot
logroot = args.logroot

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Create log directory and out directory

check_dir(logroot)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

slr_cmd = "/nfs/research2/goldman/gregs/sw/slr-1.4.1/bin/Slr {0}"

for infile in glob(path.join(ctlroot, "Eutheria", "*/*_matched.ctl")):
    print infile
    
    basename = path.basename(infile).partition('.')[0]
    basename = "".join(basename.split("_")[0] + "_" + basename.split("_")[1])
    prefix = basename.partition('_')[0][:2]    
    
    logfile = path.join(logroot, basename + '.log')   
    
    slr = slr_cmd.format(infile)
    
    if os.getcwd() != ctlroot:
        os.chdir(ctlroot)
    
    p = Popen(['bsub','-R', 'rusage[tmp=512]', '-o'+logfile, 
                   slr])
        
    p.wait()
    
        
    