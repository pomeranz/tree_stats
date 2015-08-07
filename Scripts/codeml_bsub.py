# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 10:10:55 2015

@author: gideon

Description:
This script uses the analyze paml script and applies it to every file.
"""
from subprocess import Popen
from glob import glob
import argparse
from Bio.Phylo.PAML import codeml
import os
from os import path
from utils import check_dir


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# CMD interface with argparse

argparser = argparse.ArgumentParser()

argparser.add_argument("--alignroot", metavar="Alignment Dir", type=str, required=True)
argparser.add_argument("--treeroot", metavar="Tree Dir", type=str, required=True)
argparser.add_argument("--outroot", metavar="Output Dir + working dir + logdir", type=str, required=True)

args = argparser.parse_args()

alignroot = args.alignroot
treeroot = args.treeroot
outroot = args.outroot

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Create log directory and out directory

check_dir(outroot)

out_pre_dir = path.join(outroot, "out")
check_dir(out_pre_dir)

logdir = path.join(outroot, "logs")
check_dir(logdir)

if os.getcwd() != outroot:
    os.chdir(outroot)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# define the paml command
codeml_cmd = "python /nfs/research2/goldman/pomeranz/tree_stats/Scripts/analyse_codeml.py --alignfile {0} --treefile {1} --template_dir {2} --outfile {3} --workdir {4}"
# start the loop

for infile in glob(path.join(alignroot, "*/*.phy")):
    print infile
    
    

    os.system("sed '1 s/\([[:space:]]*[0-9]*\)\([[:space:]]*[0-9]*\)\([[:space:]]*I\)*/\1\2 I/g' {0} > {0}_tmp && mv {0}_tmp {0}".format(infile))
    
    basename = path.basename(infile).partition('.')[0]
    basename = "".join(basename.split("_")[0] + "_" + basename.split("_")[1])
    prefix = basename.partition('_')[0][:2]
    
    logfile = path.join(logdir, basename + '.log')    
    
    out_sub_dir = path.join(out_pre_dir, prefix)
    check_dir(out_sub_dir)
    outdir = path.join(out_sub_dir, basename)
    check_dir(outdir)
    outfile = path.join(outdir, basename + ".mlc")
    # fixed_outfile = path.join(outdir, basename, ".fixed.mlc")
    
    treedir = path.join(treeroot, prefix)
    treefile = path.join(treedir, basename + '.nh')
    
    # template_file = "template.ctl"
    template_dir = outroot
    
    if os.getcwd() != outdir:
        os.chdir(outdir)
    
    paml = codeml_cmd.format(infile, treefile, template_dir, outfile, outdir)
    
    p = Popen(['bsub','-R', 'rusage[tmp=512]', '-o'+logfile, '-g', '/codeml', 
                   paml])
        
    p.wait()
    
    

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
