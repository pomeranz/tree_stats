# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 16:28:29 2015

@author: gideon

###                         TREE_STATS_V2                                   ###

Description:

This script takes all the trees from the directory on 
ebi-004 /nfs/research2/goldman/gregs/slr_pipeline/data/ens/78/trees
It applies a few statistics onto it including:

-total tree length
-minimum/maximum branch length
-number of leaf nodes
-number of different species
-number of paralogs 
(that is number of leaf nodes - number of different species)
-number of human sequences


It also stores them in a tabular file which can then be accessed.

Instructions: 

Replace the directory variable to apply it to any directroy.
Change FileName, to name the output file.

"""

## INPUT 

Directory = "/home/gideon/Documents/mphil_internship/Trees/*/*" 

# It will save the file in the current directory
FileName = "tree_stats_output.csv"

#-----------------------------------------------------------------------------#

## PACKAGES

# to change directories etc..
import os
# Package to access the files on the server. 
import glob
# import regular expressions
import re
# import Tree module
from ete2 import Tree
# for writing to file
import csv

#-----------------------------------------------------------------------------#

## LOOP PREPARATION

# create a regexp to match for later
ens_RE = re.compile("[A-Z]*")

# tree list to hold the final output. 
tree_list = list() 

#-----------------------------------------------------------------------------#

## LOOP

for p in glob.glob(Directory):
    
    # list for that particular tree
    tree = list()
    
    # acces the directory of the tree file
    current_tree_directory = p
    
    # create ete tree object
    current_tree = Tree(newick = current_tree_directory)
    
    # Add tree directory for identification
    tree.append(current_tree_directory)
    
    
    ## TREE LENGTH + MAX/MIN BRANCH LENGTH
    
    max_dist = 0.0
    tree_length = 0.0
    for n in current_tree.traverse():
        tree_length += n.dist
        if n.dist > max_dist:
            max_dist = n.dist
                        
    # add tree length
    tree.append(tree_length)
    
    # add max branch length
    tree.append(max_dist)
    
    # calculate min dist
    
    min_dist = 10000.0
    for n in current_tree.traverse():
        if n.dist < min_dist:
            min_dist = n.dist
    
    # add minimum branch length
    tree.append(min_dist)
    
    
    ## MAX/MIN BRANCH LENGTHS FROM ROOT
    
    # max length
    max_leaf = current_tree.get_farthest_leaf()
    #add to list
    tree.append(max_leaf)
    
    # min length
    min_leaf = current_tree.get_closest_leaf()
    # add to list
    tree.append(min_leaf)
    
    
    # NUMBER OF LEAVES
    
    # calculate number of leaves
    no_leaves = len(current_tree.get_leaves())
    # add info to tree list
    tree.append(no_leaves)
    
    
    # NUMBER OF DIFFERENT SPECIES
    
    # save all the names in an object
    leaf_names = current_tree.get_leaf_names()
    
    # use regexp to extract only species ids
    species_ids = [ ens_RE.match(s).group(0) for s in leaf_names ]
    
    unique_species_ids = list(set(species_ids))
    
    no_species = len(unique_species_ids)
    
    # add to list
    tree.append(no_species)
    
    
    ## NUMBER OF PARALOGS
    
    # paralogs are number of leaves - number of sepcies
    no_paralogs = no_leaves - no_species
    
    # add to list
    tree.append(no_paralogs)
    
    ## NUMBER OF HUMAN SEQUENCES
    
    # count the number of only ENS for human seqs
    human_seqs = species_ids.count("ENS")
    
    if human_seqs == 0:
        no_human_seqs = 0
        tree.append(no_human_seqs)
    else:
        no_human_seqs = len(human_seqs)
        tree.append(no_human_seqs)
    
    
    ## ADD TREE TO TREE LIST
    tree_list.append(tree)
    
    ## END OF LOOP
    
#-----------------------------------------------------------------------------#

## WRITE TO FILE

# create the column names
colnames = ["Directory", "Length", "Biggest Branch", "Smallest Branch", 
              "Farthest Leaf", "Closest Leaf", "# of Leaves", "# of Species", 
              "# of Paralogs", "# of Human seqs"]

# wrtite it to a file
with open(FileName, 'w') as csvfile:
    writer = csv.writer(csvfile, delimiter=",", dialect="excel")
    writer.writerow(colnames)
    writer.writerows(tree_list)
