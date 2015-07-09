# -*- coding: utf-8 -*-
"""
Created on Thu Jul  9 16:27:16 2015

@author: gideon
"""

from Bio import SeqIO



file = "/home/gideon/Documents/mphil_internship/prank_out/11_1_prank.best.fas"

fast1 = SeqIO.parse(file, "fasta")

for seq in SeqIO.parse(file, "fasta"):
        
        print seq.id
        #SeqIO.write(seq_record, "/Users/bucephalus/Desktop/nc_pimps/raw/"+name[0]+".fas", "fasta")
        
"""
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