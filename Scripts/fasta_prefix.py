# -*- coding: utf-8 -*-
"""
Created on Thu Jul  9 16:27:16 2015

@author: gideon
"""

from Bio import SeqIO
from argparse import ArgumentParser
from ete2 import Tree
# import re 
from glob import glob


argparser = ArgumentParser()

# commands for function
argparser.add_argument("--fasta", metavar="Fasta files directory", type=str, required=True)
argparser.add_argument("--tree", metavar="Tree files directory", type=str, required=True)


def fasta_prefix(fasta, tree):
    
    tree_dir = ''.join([tree + "*/*.nh"])
    fasta_dir = ''.join([fasta + "*/*prank.best.fas"])
    
    tree_files = sorted(glob(tree_dir))
    fasta_files = sorted(glob(fasta_dir))
    
    for file in range(len(fasta_files)-1):
        
        # print progress
        print("".join("Finished " + str(file + 1) + " out of " + str(len(fasta_files))))        
        
        current_fasta = fasta_files[file]
        current_tree = tree_files[file]
        
        #current_fasta = SeqIO.parse(current_fasta,"fasta")
        
        for ID in SeqIO.parse(current_fasta,"fasta"):
            # create a regexp to search for in the other file
            # id_RE = re.compile(ID)
            
            tree_ids = Tree(newick=current_tree)
            
            for tree_id in tree_ids.iter_leaf_names():
                
                if ID.id in tree_id:
                    ID.id.replace(ID.id, tree_id)
            
        SeqIO.write(current_fasta, current_fasta, "fasta")   
        
args = argparser.parse_args()

fasta = args.fasta
tree = args.tree

fasta_prefix(fasta,tree)

print("All done, the fasta file names should match the tree IDs")




"""
for ID in SeqIO.parse(current_fasta,"fasta"):
     tree_ids = Tree(newick=current_tree)
     for tree_id in tree_ids.iter_leaf_names():
         print ID.id
         if ID.id in tree_id:
             ID.id = tree_id
             print(ID.id) 






file = "/home/gideon/Documents/mphil_internship/prank_out/11_1_prank.best.fas"

fast1 = SeqIO.parse(file, "fasta")

for seq in SeqIO.parse(file, "fasta"):
        
        print seq.id
        #SeqIO.write(seq_record, "/Users/bucephalus/Desktop/nc_pimps/raw/"+name[0]+".fas", "fasta")
       

Read in both directories
Within directory read in correct files. 1_1.nh and 1_1.fas

create a list with all the leaf names from the nh file.

Loop trough fasta headers.

for each header create a regexp and search for it within the list.
The regexp greps out the 3 characters before the match = D1_ or D2_ 
add this to the current header.
Next header.

Next fasta and nh file
Next directory.

Done!
"""