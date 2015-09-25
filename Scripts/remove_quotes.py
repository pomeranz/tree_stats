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
    
    tree_file = open(infile, "r+")
    
    text = tree_file.read()
    text = text.replace("'", "")
    tree_file.seek(0)
    tree_file.write(text)
    tree_file.truncate()
    tree_file.close()
    
    tree_file = path.join(args.preproot, "*", "*.nwk")   
    
    tree = Tree.get_from_path(tree_file, 'newick', preserve_underscores=True)
    
    for node in tree:
        if "." in node.taxon.label:
            node.taxon.label = node.taxon.label.replace(".", "_")
    
    tree_file.seek(0)
    tree_file.write(tree)
    tree_file.truncate()
    tree_file.close
    
    
    fasta_file = path.join(args.preproot, "*", "*.fa")
    
    text = fasta_file.read()
    text = text.replace(".", "")
    fasta_file.seek(0)
    fasta_file.write(text)
    fasta_file.trunctuate()
    fasta_file.close
