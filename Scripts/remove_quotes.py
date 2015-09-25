# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 16:31:04 2015

@author: gideon

Description: This script removes all quotes from the raxml trees.
"""

from glob import glob
import os
from os import path
import argparse
from dendropy import Tree


argparser = argparse.ArgumentParser()

argparser.add_argument('--preproot', metavar='tree_root', type=str, required=True)


args = argparser.parse_args()

for infile in glob(path.join(args.preproot, "*", "*.nwk")):
    
    print infile
    
    basename = path.basename(infile).partition('.')[0]
    prefix = basename.partition('_')[0][:2]
    
    tree_file = open(infile, "r+")
    
    text = tree_file.read()
    text = text.replace("'", "")
    tree_file.seek(0)
    tree_file.write(text)
    tree_file.truncate()
    tree_file.close()
       
    
    tree = Tree.get_from_path(infile, 'newick', preserve_underscores=True)
    
    for node in tree:
        if node.is_leaf():
            if "." in node.taxon.label:
                node.taxon.label = node.taxon.label.replace(".", "_")
            
    tree_file = open(infile, "r+")
    
    tree_file.seek(0)
    tree_file.write(tree.as_string('newick'))
    tree_file.truncate()
    tree_file.close
    
    
    fasta_file = open(path.join(args.preproot, prefix, basename + ".fa"), "r+")
    
    text = fasta_file.read()
    text = text.replace(".", "")
    fasta_file.seek(0)
    fasta_file.write(text)
    fasta_file.trunctuate()
    fasta_file.close
