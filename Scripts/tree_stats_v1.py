# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 14:37:36 2015

@author: gideon

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

It should also save all the trees so you can have a look at all of them,
if one wants
"""
# to change directories etc..
import os

# import ete2 to visualize trees
from ete2 import EvolTree
# import special representation of ete2
from ete2.treeview.layouts import evol_clean_layout

# Package to access the files on the server. 
import glob

# Another tree manipulation package
import dendropy

# import regular expressions
import re

"""
Snippet from Greg that lets you print out the directory names of the trees

import glob
 
for d in glob.glob("/nfs/research2/goldman/gregs/slr_pipeline/data/ens/78/trees/*"):
  print "Directory name", d
 
 
# Probably more useful but would print a lot of stuff
# for f in glob.glob("/nfs/research2/goldman/gregs/slr_pipeline/data/ens/78/trees/*/*.nh"):
#   print "File name", f
 
"""

# import Tree module
from ete2 import Tree

# test a tree

tree_1 = Tree("/home/gideon/Documents/mphil_internship/Trees/1/1.nh")
print tree_1
tree_2 = Tree("/home/gideon/Documents/mphil_internship/Trees/2/2.nh")
print tree_2

# test how to show a tree
from ete2 import Tree, TreeStyle

t = Tree( "((a,b),c);" )
circular_style = TreeStyle()
circular_style.mode = "c" # draw tree in circular mode
circular_style.scale = 20
t.render("mytree.png", w=183, units="mm", tree_style=circular_style)

t = Tree("/home/gideon/Documents/mphil_internship/Trees/1/1.nh")
t.populate(10, random_dist=True)
ts = TreeStyle()
ts.show_leaf_name = True
ts.show_branch_length = True
ts.show_branch_support = True
t.show(tree_style=ts)




#------------------------------------------------------------------#

### FUNCTIONS ###


## Tree length

# create a tree object
tree_dendro = dendropy.Tree.get(
path="/home/gideon/Documents/mphil_internship/Trees/1/1.nh", 
schema="newick")  # specify path + schema

# get length of tree
tree_dendro.length()


## Max/Min branch from root

# use ete2 object
from ete2 import Tree

# create obejct
tree_ete = Tree(newick="/home/gideon/Documents/mphil_internship/Trees/1/1.nh")

# get closest leaf
tree_ete.get_closest_leaf()

# get furthest tree
tree_ete.get_farthest_leaf()

## Max/Min branch length in the tree + tree length
max_dist = 0.0
tree_length = 0.0
for n in tree_ete.traverse():
    tree_length += n.dist
    if n.dist > max_dist:
        max_dist = n.dist



## Number of Leaf nodes
from ete2 import Tree
tree_ete = Tree(newick="/home/gideon/Documents/mphil_internship/Trees/1/1.nh")

no_leaves = len(tree_ete.get_leaves())

## Number of species

# get names of all the leaves
tree_ete.get_leaf_names()

# save all the names in an object
leaf_names = tree_ete.get_leaf_names()

# keep only the species identifier
leaf_names2 = [leaf_names[x][:7] for x in range(0,len(leaf_names))]

# convert it to a set which gives unique names and back into a list
leaf_names3 = list(set(leaf_names2))

# get number of entries = number of species
len(leaf_names3)

## Number of human sequences





#------------------------------------------------------------------#

"""

# This was the first try. Actual loop further down. 

# This is a loop which is going to create the tree objects.
# Later to be included in big loop


# tree list to hold the final output. 
tree_list = list()


for p in glob.glob("/home/gideon/Documents/mphil_internship/Trees/*/*"):
    
    # acces the directory of the tree file
    current_tree_directory = p
    # create the tree object
    current_tree = dendropy.Tree.get(path=current_tree_directory, 
                                     schema="newick")
    # add tree to list for later
    tree_list.append(current_tree)

"""

# tree list to hold the final output. 
tree_list = list() 

# create a regexp to match for later
ens_RE = re.compile("[A-Z]*")
    
#  Test loop that iterates over the trees and collects the statistics

# Each tree is going to be a list, with its stats. Then we add the tree 
# to a list of trees.
# In the end we can just ask what is the longest tree, etc...


for p in glob.glob("/home/gideon/Documents/mphil_internship/Trees/*/*"):
    
    # list for that particular tree
    tree = list()
    
    # acces the directory of the tree file
    current_tree_directory = p
    
    # create ete tree object
    current_tree = Tree(newick = current_tree_directory)
    
    # Add tree name
    tree.append(current_tree_directory)
    
    """
    # create the tree object
    current_tree = Tree(newick=current_tree_directory)
    
    # first entry is the tree name
    tree.append(current_tree)
    
    # TREE LENGTH
    tree_length = current_tree.length()
    # add length as second entry to the tree list
    tree.append(tree_length)
    
    """
    
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
    
    """
    # keep only the species identifier
    leaf_names2 = [leaf_names[x][:7] for x in range(0,len(leaf_names))]

    # convert it to a set which gives unique names and back into a list
    leaf_names3 = list(set(leaf_names2))

    # get number of entries = number of species
    no_species = len(leaf_names3)
    
    """
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
    
    """
    
    # That was experimenting
    
    [s for s in leaf_names2 if (["ENS" "\d*"] in s)]
    
    m = re.search('[A-Z]+', bla[0])
    m = m.group(0)
    m = m[:-1]
    
    no_human_seqs = 0
    
    for m in range(0,len(leaf_names)):
        s = leaf_names[m]
        r = re.search('ENS[0-9]', s)
        if r is not None:
            no_human_seqs = no_human_seqs + 1
            
    """
    
    ## ADD TREE TO TREE LIST
    tree_list.append(tree)
    
    ## END OF LOOP
    
    
#-------------------------------------------------------------#
    
# This is about writing everything to a tabular file

import csv

colnames = ["Directory", "Length", "Biggest Branch", "Smallest Branch", 
              "Farthest Leaf", "Closest Leaf", "# of Leaves", "# of Species", 
              "# of Paralogs", "# of Human seqs"]


with open('test_file.csv', 'w') as csvfile:
    writer = csv.writer(csvfile, delimiter=",", dialect="excel")
    writer.writerow(colnames)
    writer.writerows(tree_list)
