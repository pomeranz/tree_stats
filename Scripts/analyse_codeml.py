# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 13:08:06 2015

@author: gideon

Description:
This script uses Biopythons implementation of codeml.
"""

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
    
    cml.alignment = alignfile
    #fixed_cml.alignment = alignfile
    cml.tree = treefile
    #fixed_cml.tree = treefile
    cml.out_file = "".join(outfile + ".mlc")
    #fixed_cml.out_file = "".join(outfile, ".fixed.mlc")
    cml.working_dir = workdir 
    #fixed_cml.working_dir = workdir
    
    
    cml.read_ctl_file(path.join(template_dir, "template.ctl"))
    #fixed_cml.read_ctl_file(path.join(template_dir, ".fixed.ctl"))
    
    cml.run()
    #fixed_cml.run()
    

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# call the function
analyse_codeml(alignfile, treefile, template_dir, outfile, workdir)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%