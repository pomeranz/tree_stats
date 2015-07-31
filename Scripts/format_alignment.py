# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 11:21:04 2015

@author: gideon

Description:
This script takes the output of PRANK and formats everything
so that the alignment can be used for tdg09
- Prefixes
- Format conversion (fasta -> phylip-relaxed)
- 4 alignments available
    - fasta
    - fasta AA
    - phylip
    - phylip AA

Usage:
python Scripts/format_align.py --fasta --tree --outdir


To test:
tree = "/home/gideon/Documents/mphil_internship/fas_pref_test/Eutheria/"
fasta = "/home/gideon/Documents/mphil_internship/fas_pref_test2/"
outdir = "/home/gideon/Documents/mphil_internship/alignments2"


"""
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
import os
from os import path
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

from utils import check_dir

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
argparser = ArgumentParser()

# commands for function
argparser.add_argument("--fasta", metavar="Fasta files directory", type=str, required=True)
argparser.add_argument("--tree", metavar="Tree files directory", type=str, required=True)
argparser.add_argument("--outdir", metavar="Directory to write the new file to", type=str, required=True)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def format_alignment(fasta, tree, outdir):
    
    treeroot = tree
    fastaroot = path.join(fasta, "*/*prank.best.fas")
    
    check_dir(outdir)


    for infile in glob(fastaroot):
        
        # print progress
        print infile
        
        basename = path.basename(infile).partition('.')[0]
        basename = "".join(basename.split("_")[0] + "_" + basename.split("_")[1])
        prefix = basename.partition('_')[0][:2]
        
        fastafile = infile 
        treedir = path.join(treeroot, prefix)
        treefile = path.join(treedir, basename + '.nh')
        
        # create the first 2 directories (fasta_out, fasta_AA_out)
        
        fasta_out_dir = path.join(outdir, "fasta")
        check_dir(fasta_out_dir)
        fasta_AA_out_dir = path.join(outdir, "fasta_AA")
        check_dir(fasta_AA_out_dir)
        
        fasta_out_subdir = path.join(fasta_out_dir, prefix)
        check_dir(fasta_out_subdir)
        fasta_out_file_path = path.join(fasta_out_subdir, "".join(basename + ".fa"))
        fasta_AA_out_subdir = path.join(fasta_AA_out_dir, prefix)
        check_dir(fasta_AA_out_subdir)
        fasta_AA_out_file_path = path.join(fasta_AA_out_subdir, "".join(basename + ".fa"))
        
        fasta_out_file = open(fasta_out_file_path, "w")
        fasta_AA_out_file = open(fasta_AA_out_file_path, "w")        

          
        for ID in SeqIO.parse(fastafile,"fasta", alphabet=IUPAC.unambiguous_dna):
            
            tree_ids = Tree(newick=treefile)
            for tree_id in tree_ids.iter_leaf_names():
                
                if tree_id.find(ID.id) != -1:
                    #print ID.id
                    ID.id = tree_id
                    #ID.name = ""
                    ID.description = ""
                    #print ID.id
                    #print ID
                    
                    # write the normal fasta out
                    SeqIO.write(ID, fasta_out_file, "fasta")
                    
                    # translate cDNA and write AA fasta
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

                    ID = SeqRecord(Seq(aa_seq, IUPAC.protein), id = ID.id, name = ID.name)
                    ID.description = ""

                    SeqIO.write(ID, fasta_AA_out_file, "fasta")
                    
        fasta_out_file.close()
        fasta_AA_out_file.close()
        
        phy_out_dir = path.join(outdir, "phylip")
        check_dir(phy_out_dir)
        phy_AA_out_dir = path.join(outdir, "phylip_AA")
        check_dir(phy_AA_out_dir)
        
        phy_out_subdir = path.join(phy_out_dir, prefix)
        check_dir(phy_out_subdir)
        phy_out_file_path = path.join(phy_out_subdir, "".join(basename + ".phy"))
        phy_AA_out_subdir = path.join(phy_AA_out_dir, prefix)
        check_dir(phy_AA_out_subdir)
        phy_AA_out_file_path = path.join(phy_AA_out_subdir, "".join(basename + ".phy"))

        fasta_alignment = open(fasta_out_file_path, "rU")
        fasta_AA_alignment = open(fasta_AA_out_file_path, "rU")
        
        phy_out_file = open(phy_out_file_path, "w")
        phy_AA_out_file = open(phy_AA_out_file_path, "w")
                        
        alignments = AlignIO.parse(fasta_alignment, "fasta")
        AlignIO.write(alignments, phy_out_file, "phylip-relaxed")
        
        fasta_alignment.close()
        phy_out_file.close()

        alignments_AA = AlignIO.parse(fasta_AA_alignment, "fasta")       
        AlignIO.write(alignments_AA, phy_AA_out_file, "phylip-relaxed")

        fasta_AA_alignment.close()
        phy_AA_out_file.close()
           
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
args = argparser.parse_args()

fasta = args.fasta
tree = args.tree
outdir = args.outdir


format_alignment(fasta, tree, outdir)

print("All done, the fasta file names should match the tree IDs")
print("The alignments are available in fasta, fasta(AA), phylip, phylip(AA)")
