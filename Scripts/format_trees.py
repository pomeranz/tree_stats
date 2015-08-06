# -*- coding: utf-8 -*-
"""
Created on Sun Jul 19 14:41:11 2015

@author: gideon

Description:
This script takes the tree directory and the make_seqsets_out as inputs.
It iterates over the fasta files and removes tree_leaves which are NOT 
in the Fasta file.
Then it outputs the tree files to a new directory


Test:
treeroot = "/home/gideon/Documents/mphil_internship/fas_pref_test/Eutheria/"
fastaroot = "/home/gideon/Documents/mphil_internship/fas_pref_test2/"
outroot = "/home/gideon/Documents/mphil_internship/format_tree_out"

"""
import os
from os import path
from glob import glob
from ete2 import Tree
from Bio import SeqIO
from argparse import ArgumentParser
from utils import check_dir
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

argparser = ArgumentParser()

argparser.add_argument("--treeroot", metavar="Tree Directory", type=str, required=True)
argparser.add_argument("--fastaroot", metavar="Fasta Directory", type=str, required=True)
argparser.add_argument("--outroot", metavar="Output Directory", type=str, required=True)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def format_trees(treeroot, fastaroot, outroot):
    
    fastafiles = path.join(fastaroot, "*/*.fa")
    
    if not os.path.exists(outroot):
        os.makedirs(outroot)
    
    rooted_out_dir = path.join(outroot, "rooted")
    check_dir(rooted_out_dir)
    unrooted_out_dir = path.join(outroot, "unrooted")
    check_dir(unrooted_out_dir)
    
    
    for infile in glob(fastafiles):
        
        print infile
        
        basename = path.basename(infile).partition('.')[0]
        basename = "".join(basename.split("_")[0] + "_" + basename.split("_")[1])
        prefix = basename.partition('_')[0][:2]
        
        fastafile = infile 
        treedir = path.join(treeroot, prefix)
        treefile = path.join(treedir, basename + '.nh')
        
        # make the tree object
        tree = Tree(newick=treefile)
        
        # loop that deletes nodes that are not in the alignment
        for leaf_name in tree.iter_leaf_names():
            
            name_check = []
            
            for ID in SeqIO.parse(fastafile, "fasta"):
                if ID.id in leaf_name:
                    name_check.append(True)
                else:
                    name_check.append(False)
            
            if any(name_check):
                continue
            else:
                leaf = tree.search_nodes(name=leaf_name)[0]
                leaf.delete()
                #node = leaf.up
                #node.remove_child(leaf)
                    
            # create the directories for rooted trees
            rooted_out_sub_dir = path.join(rooted_out_dir, prefix)
            check_dir(rooted_out_sub_dir)
            rooted_out_file = path.join(rooted_out_sub_dir, basename + ".nh")
            
            
            
            tree.write(outfile=rooted_out_file, format=6)
            
            # create subdirectories for unrooted trees
            unrooted_out_sub_dir = path.join(unrooted_out_dir, prefix)
            check_dir(unrooted_out_sub_dir)
            unrooted_out_file = path.join(unrooted_out_sub_dir, basename + ".nh")
            # unroot the tree
            tree.unroot()
            
            tree.write(outfile=unrooted_out_file, format=6)
            
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# call the function
args = argparser.parse_args()

treeroot = args.treeroot
fastaroot = args.fastaroot
outroot = args.outroot

format_trees(treeroot, fastaroot, outroot)

print("All done, the trees should be formatted now. 2 directories: rooted and unrooted created!")