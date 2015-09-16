# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 10:52:24 2015

@author: gideon
"""
import time
from subprocess import Popen
from glob import glob
import argparse
from Bio.Phylo.PAML import codeml
import os
from os import path
from utils import check_dir
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

if os.getcwd() != out_pre_dir:
    os.chdir(out_pre_dir)
    

for infile in glob(path.join(alignroot, "*/*.phy")):
    print infile
    
    with open(infile, 'r') as f:
        first_line = f.readline()
        if "I" not in first_line:
            os.system("sed '1 s/$/ I/' {0} > {0}_tmp && mv {0}_tmp {0}".format(infile))
    
    basename = path.basename(infile).partition('.')[0]
    basename = "".join(basename.split("_")[0] + "_" + basename.split("_")[1])
    prefix = basename.partition('_')[0][:2]
    
    treedir = path.join(treeroot, prefix)
    treefile = path.join(treedir, basename + '.nh')
    
    
                    #### OLD PAML RUN ###
    old_human_time_start = time.localtime()    
    
    old_start_time = time.time()    
    
    # out directories
    old_out_pre_sub_dir = path.join(out_pre_dir,"old")
    check_dir(old_out_pre_sub_dir)
    old_out_sub_dir = path.join(old_out_pre_sub_dir, prefix)
    check_dir(old_out_sub_dir)   
    
    old_outdir = path.join(old_out_sub_dir, basename)
    check_dir(old_outdir)
    old_outfile = path.join(old_outdir, basename + ".mlc")
    
    if os.getcwd() != old_outdir:
        os.chdir(old_outdir)
    
    old_cml = codeml.Codeml()
    old_cml.read_ctl_file("template.ctl")
    old_cml.alignment = infile
    old_cml.tree = treefile
    old_cml.out_file = old_outfile
    old_cml.working_dir = old_outdir
    
    # create ctl file
    old_ctlfile = "".join(old_outfile.split(".")[0] + ".ctl")
    old_cml.ctl_file = old_ctlfile
    old_cml.write_ctl_file()
    old_cml.run(ctl_file=old_ctlfile, command = "/nfs/research2/goldman/botond/soft/bin/codeml", verbose=True)
    
    old_end_time = time.time()
    old_human_time_end = time.localtime()
    # runtime in seconds
    old_runtime = old_end_time - old_start_time
    
    # write time to a file
    
    old_times = open(path.join(old_outdir,"times.txt"), "w")
    old_times.write(old_human_time_start, old_human_time_end, old_start_time, old_end_time, old_runtime)
    
    

                    #### NEW PAML RUN ###
    
    new_human_time_start = time.localtime()    
    
    new_start_time = time.time()
    
    new_out_pre_sub_dir = path.join(out_pre_dir,"new")
    check_dir(new_out_pre_sub_dir)
    new_out_sub_dir = path.join(new_out_pre_sub_dir, prefix)
    check_dir(new_out_sub_dir)
    new_outdir = path.join(new_out_sub_dir, basename)
    check_dir(new_outdir)
    new_outfile = path.join(new_outdir, basename + ".mlc")
    
    if os.getcwd() != new_outdir:
        os.chdir(new_outdir)
    
    new_cml = codeml.Codeml()
    new_cml.read_ctl_file("template.ctl")
    new_cml.alignment = infile
    new_cml.tree = treefile
    new_cml.out_file = old_outfile
    new_cml.working_dir = old_outdir
    

    # create ctl file
    new_ctlfile = "".join(new_outfile.split(".")[0] + ".ctl")
    new_cml.ctl_file = new_ctlfile
    new_cml.write_ctl_file()
    new_cml.run(ctl_file=new_ctlfile, command = "/nfs/research2/goldman/gregs/sw/paml4.8/bin/codeml", verbose=True)
    
    new_end_time = time.time()
    new_human_time_end = time.localtime()
    # time in seconds
    new_runtime = new_end_time - new_start_time
    
    # write time to a file # finish this
    new_times = open(path.join(new_outdir,"times.txt"), "w")
    new_times.write(new_human_time_start, new_human_time_end, new_start_time, new_end_time, new_runtime)
    
    