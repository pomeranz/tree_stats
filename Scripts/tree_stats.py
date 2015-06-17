# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 16:49:53 2015

@author: gideon
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 16:28:29 2015

@author: gideon

###                         TREE_STATS                                   ###


Description:

This script defines a function that collects some statistics from a tree 
directory. The directory needs to be passed to the function. 

-total tree length
-minimum/maximum branch length
-number of leaf nodes
-number of different species
-number of paralogs 
(that is number of leaf nodes - number of different species)
-number of human sequences


It also stores them in a tabular file which can then be accessed.

Instructions: 

Run the script. It will create a function.
It needs to be passed a directory with the locations of 
the tree files.
Optional arguments are the ouput file name and the directory
in which the results should be put. 

"""

## INPUT 

# directory = "/home/gideon/Documents/mphil_internship/Trees/*/*" 

# It will save the file in the current directory
# filename = "tree_stats_output.csv"
    

    
    
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
# for creating command line interface
import argparse

#-----------------------------------------------------------------------------#

### THE FUNCTION ###

def tree_stats(directory, filename="tree_stats_output.csv", 
               output_directory = os.getcwd()):
    
    ## LOOP PREPARATION
    
    # match all the sub-directories in the directory
    directory = ''.join([directory + "/*/*"])

    # create a regexp to match for later
    ens_RE = re.compile("[A-Z]*")

    # tree list to hold the final output. 
    tree_list = list() 
    
#-----------------------------------------------------------------------------#
    
    ## LOOP
    
    for p in glob.glob(directory):
        
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
        
        ## Show progress
        print("Current file:" + current_tree_directory)        
        
        
        ## END OF LOOP
        
#-----------------------------------------------------------------------------#
    
    ## WRITE TO FILE
    
    # create the column names
    colnames = ["Directory", "Length", "Biggest Branch", "Smallest Branch", 
                  "Farthest Leaf", "Closest Leaf", "# of Leaves", 
                  "# of Species", "# of Paralogs", "# of Human seqs"]
                  
    output_file = ''.join([output_directory + "/" + filename ])
    
    # wrtite it to a file
    with open(output_file, 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter=",", dialect="excel")
        writer.writerow(colnames)
        writer.writerows(tree_list)
            
    ## END OF FUNCTION
    
#-----------------------------------------------------------------------------#

# create the argparse so that it can be run from the command-line. 
    
parser = argparse.ArgumentParser(prog="tree_stats",
                                 description='Calculates tree statistics')

# add directroy argument, must be given
parser.add_argument("--directory", metavar="Tree Directory", type=str)

# give output file name, optional
parser.add_argument("--outname", metavar="Output file name", action="store", 
                    default="tree_stats_output.csv", required=False,
                    dest="filename")

# add optional argument for output directory
parser.add_argument("--outdir", metavar="Output Directory", action="store",
                    default=os.getcwd(), required=False, 
                    dest="output_directory")

args = parser.parse_args()


## CALL THE FUNCTION

tree_stats(args.directory)

## TELL THE USER WE ARE DONE
print("Finished, the file output is in your current/specified working directory")
