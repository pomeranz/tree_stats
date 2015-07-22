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
from glob import glob
from ete2 import Tree
from Bio import SeqIO
from argparse import ArgumentParser

argparser = ArgumentParser()

argparser.add_argument("--treeroot", metavar="Tree Directory", type=str, required=True)
argparser.add_argument("--fastaroot", metavar="Fasta Directory", type=str, required=True)
argparser.add_argument("--outroot", metavar="Output Directory", type=str, required=True)


def format_trees(treeroot, fastaroot, outroot):
    
    tree_dir = ''.join([treeroot + "*/*.nh"])
    fasta_dir = ''.join([fastaroot + "*/*.fa"])
    
    if not os.path.exists(outroot):
        os.makedirs(outroot)
    
    tree_files = sorted(glob(tree_dir))
    fasta_files = sorted(glob(fasta_dir))
     
    
    for index in range(len(fasta_files)):
        current_fasta = fasta_files[index]
        current_tree = tree_files[index]
        
        # this is wrong! need to only remove the tree leaves
        # which are not in the fasta file  
        tree = Tree(newick=current_tree)
        
        for leaf_name in tree.iter_leaf_names():
            
            name_check = []
            
            for ID in SeqIO.parse(current_fasta, "fasta"):
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
                    
            
            out_sub_dir = "".join(outroot + "/" + current_tree.split("/")[7])
            out_file = "".join(out_sub_dir + "/" + current_tree.split("/")[8])
            
            if not os.path.exists(out_sub_dir):
                os.makedirs(out_sub_dir)
            
            
            tree.write(outfile=out_file, format=6)
        print("".join("Finished " + str(index + 1) + " out of " + str(len(fasta_files))))        
            

# call the function
args = argparser.parse_args()

treeroot = args.treeroot
fastaroot = args.fastaroot
outroot = args.outroot

format_trees(treeroot, fastaroot, outroot)

print("All done, the trees should be formatted now")