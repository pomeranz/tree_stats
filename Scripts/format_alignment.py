# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 11:21:04 2015

@author: gideon

Description:
This script takes the output of PRANK and formats everything
so that the alignment can be used for tdg09
- Prefixes
- Format conversion (fasta -> phylip-relaxed)

Side note: It creates an intermediate directory with just the prefixes 
changed. 

Usage:
python Scripts/format_align.py --fasta --tree --outdir --final_out

To test:
tree = "/home/gideon/Documents/mphil_internship/fas_pref_test/Eutheria/"
fasta = "/home/gideon/Documents/mphil_internship/fas_pref_test2/"
outdir = "/home/gideon/Documents/mphil_internship/fasta_prefix_out"
final_out = "/home/gideon/Documents/mphil_internship/alignments2"
toAA = True
"""

import os
from Bio import SeqIO
from argparse import ArgumentParser
from ete2 import Tree
# import re
from glob import glob
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Seq import translate
import itertools
from itertools import izip_longest
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord

def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)


argparser = ArgumentParser()

# commands for function
argparser.add_argument("--fasta", metavar="Fasta files directory", type=str, required=True)
argparser.add_argument("--tree", metavar="Tree files directory", type=str, required=True)
argparser.add_argument("--outdir", metavar="Directory to write the new file to", type=str, required=True)
argparser.add_argument("--final_out", metavar="Directory to write the new file to", type=str, required=True)
argparser.add_argument("--translate", metavar="Translate DNA?", type=str, required=True)


def fasta_prefix(fasta, tree, outdir, final_out, toAA):
    
    tree_dir = ''.join([tree + "*/*.nh"])
    fasta_dir = ''.join([fasta + "*/*prank.best.fas"])
    
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    
    tree_files = sorted(glob(tree_dir))
    fasta_files = sorted(glob(fasta_dir))
    
    # make sure that both directories have the same amount of files
    if len(tree_files) != len(fasta_files):
        index = 0
        while index <= len(tree_files)-1:
            current_fasta = fasta_files[index]
            current_fasta = current_fasta.split("/")
            current_fasta = current_fasta[7]
            current_fasta = current_fasta.split("_prank")
            current_fasta = current_fasta[0]
            
            current_tree = tree_files[index]
            current_tree = current_tree.split("/")
            current_tree = current_tree[len(current_tree)-1]
            current_tree = current_tree.split(".")
            current_tree = current_tree[0]
            if current_fasta != current_tree:
                tree_files.remove(tree_files[index])
                index = index - 1
            print index
            index += 1
    
    #tree_files = sorted(glob(tree_dir))
    fasta_files = sorted(glob(fasta_dir))
    
    for file in range(len(fasta_files)):
        
        # print progress
        print("".join("Finished " + str(file + 1) + " out of " + str(len(fasta_files))))        
        
        current_fasta = fasta_files[file]
        current_tree = tree_files[file]
        
        out_sub_dir = "".join(outdir + "/" + current_fasta.split("/")[6])
        
        if not os.path.exists(out_sub_dir):
            os.makedirs(out_sub_dir)
        
        out_file_temp = "".join(out_sub_dir + "/" + current_fasta.split("/")[7])
        out_file = open(out_file_temp, "w+")
                  
        for ID in SeqIO.parse(current_fasta,"fasta", alphabet=IUPAC.unambiguous_dna):
            
            tree_ids = Tree(newick=current_tree)
            for tree_id in tree_ids.iter_leaf_names():
                
                if tree_id.find(ID.id) != -1:
                    print ID.id
                    ID.id = tree_id
                    #ID.name = ""
                    ID.description = ""
                    print ID.id
                    print ID
                    if toAA == "False":                    
                        SeqIO.write(ID, out_file, "fasta")
                    else:
                        aa_seq = []
                        coding_dna = ID.seq
                        #print coding_dna
                        for codon in grouper(coding_dna, 3):
                            cog = "".join(codon)
                            if cog == "---":
                                aa_seq.append("-")
                            else:
                                cog_aa = translate(cog)
                                aa_seq.append(cog_aa)
                        aa_seq = "".join(aa_seq)
                        #print aa_seq
                        
                        #ID.seq = "".join(aa_seq + ", IUPACProtein()")
                        ID = SeqRecord(Seq(aa_seq, IUPAC.protein), id = ID.id, name = ID.name)
                        print ID
                        SeqIO.write(ID, out_file, "fasta")
                        
            
                        
            
        #fas.close()
        out_file.close()
        
        input_handle = open(out_file_temp, "rU")
        out_sub_dir = "".join(final_out + "/" + current_fasta.split("/")[6])
        
        if not os.path.exists(out_sub_dir):
            os.makedirs(out_sub_dir)
        
        out_file = "".join(out_sub_dir + "/" + current_fasta.split("/")[7])
        output_handle = open(out_file, "w+")
        
        alignments = AlignIO.parse(input_handle, "fasta")
        AlignIO.write(alignments, output_handle, "phylip-relaxed")
 
        output_handle.close()
        input_handle.close()
    
        
args = argparser.parse_args()

fasta = args.fasta
tree = args.tree
outdir = args.outdir
final_out = args.final_out
toAA = args.translate

fasta_prefix(fasta, tree, outdir, final_out, toAA)

print("All done, the fasta file names should match the tree IDs")
