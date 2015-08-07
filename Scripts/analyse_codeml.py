# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 13:08:06 2015

@author: gideon

Description:
This script uses Biopythons implementation of codeml.
"""
import os
from Bio.Phylo.PAML import codeml
from os import path
import argparse

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
argparser = argparse.ArgumentParser()

argparser.add_argument("--alignfile", metavar="Alignment Dir", type=str, required=True)
argparser.add_argument("--treefile", metavar="Tree Dir", type=str, required=True)
argparser.add_argument("--outfile", metavar="Output Dir + working dir + logdir", type=str, required=True)
argparser.add_argument("--template_dir", metavar="template Dir", type=str, required=True)
argparser.add_argument("--workdir", metavar="workdir Dir", type=str, required=True)

args = argparser.parse_args()

alignfile = args.alignfile
treefile = args.treefile
outfile = args.outfile
template_dir = args.template_dir
workdir = args.workdir
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def analyse_codeml(alignfile, treefile, template_dir, outfile, workdir):
    
    cml = codeml.Codeml()
    #fixed_cml = codeml.Codeml()
    
    if os.getcwd() != workdir:
        os.chdir(workdir)    
    
    cml.read_ctl_file("template.ctl")
    #fixed_cml.read_ctl_file(path.join(template_dir, ".fixed.ctl"))
    
    cml.alignment = alignfile
    print cml.alignment
    #fixed_cml.alignment = alignfile
    cml.tree = treefile
    print cml.tree
    #fixed_cml.tree = treefile
    cml.out_file = outfile
    print cml.out_file
    #fixed_cml.out_file = "".join(outfile, ".fixed.mlc")
    cml.working_dir = workdir
    print cml.working_dir
    #fixed_cml.working_dir = workdir
    
    if os.getcwd() != workdir:
        os.chdir(workdir)
    
    ctlfile = "".join(outfile.split(".")[0] + ".ctl")
    
    cml.ctl_file = ctlfile
    print cml.ctl_file
    cml.write_ctl_file()
    
    cml.run(ctl_file=ctlfile, command = "/nfs/research2/goldman/botond/soft/bin/codeml")
    #fixed_cml.run()
    

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# call the function
analyse_codeml(alignfile, treefile, template_dir, outfile, workdir)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%