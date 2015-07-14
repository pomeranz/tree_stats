# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 11:21:04 2015

@author: gideon

Rewrite:
Add another file where you write the output too. 


"""

import os
from Bio import SeqIO
from argparse import ArgumentParser
from ete2 import Tree
# import re
from glob import glob


argparser = ArgumentParser()

# commands for function
argparser.add_argument("--fasta", metavar="Fasta files directory", type=str, required=True)
argparser.add_argument("--tree", metavar="Tree files directory", type=str, required=True)
argparser.add_argument("--outdir", metavar="Directory to write the new file to", type=str, required=True)


def fasta_prefix(fasta, tree, outdir):
    
    tree_dir = ''.join([tree + "*/*.nh"])
    fasta_dir = ''.join([fasta + "*/*prank.best.fas"])
    
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    
    tree_files = sorted(glob(tree_dir))
    fasta_files = sorted(glob(fasta_dir))
    
    for file in range(len(fasta_files)-1):
        
        # print progress
        print("".join("Finished " + str(file + 1) + " out of " + str(len(fasta_files))))        
        
        current_fasta = fasta_files[file]
        current_tree = tree_files[file]
        
        out_sub_dir = "".join(outdir + current_fasta.split("/")[6])
        
        if not os.path.exists(out_sub_dir):
            os.makedirs(out_sub_dir)
        
        out_file = "".join(out_sub_dir + "/" + current_fasta.split("/")[7])
        out_file = open(out_file, "w+")
        
        with open (current_fasta, "r") as fas:
            for line in fas:
                if line.startswith(">"):
                    line = line[1:].rstrip() # skip > character
                    
                    tree_ids = Tree(newick=current_tree)
                    
                    for tree_id in tree_ids.iter_leaf_names():
                        #print line
                        #print tree_id
                        
                        
                        if tree_id.find(line) != -1:
                        # if line in tree_id:
                            print line
                            print tree_id                            
                            print True
                            new_line = tree_id
                            
                            out_file.write("".join(">" + new_line + "\n"))
                            #print line
                            #line.replace(line, new_line)
                            #break
                else:
                    out_file.write(line)
        fas.close()
        out_file.close()
        
args = argparser.parse_args()

fasta = args.fasta
tree = args.tree
outdir = args.outdir

fasta_prefix(fasta,tree, outdir)

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