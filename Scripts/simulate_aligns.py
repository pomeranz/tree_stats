# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 16:40:09 2015

@author: gideon



# For testing:
treedir = "/home/gideon/Documents/mphil_internship/simulations/test_trees/"
outdir = "/home/gideon/Documents/mphil_internship/simulations/aln1/"
gpf = "/home/gideon/Documents/mphil_internship/simulations/run1.gpf"

"""

# Packages
import pyvolve
import os
from os import path
from utils import check_dir
import argparse
import re

argparser = argparse.ArgumentParser()

argparser.add_argument("--treedir", metavar="Tree Dir", type=str, required=True)
argparser.add_argument("--outdir", metavar="Output Dir + working dir + logdir", type=str, required=True)
argparser.add_argument("--gpf", metavar="template Dir", type=str, required=True)

args = argparser.parse_args()

treedir = args.treedir
outdir = args.outdir
gpf = args.gpf


####### Prep ###########

sizes = "Small", "Medium", "Big"
species_numbers = "6species", "12species", "17species", "44species"

check_dir(outdir)
os.chdir(outdir)

# extract information out of the gpf (gideon pomeranz file)
parameters = open(gpf).read()
m = re.search("(?<=n_sites=)\w+", parameters)
n_sites = m.group(0)
n_sites = int(n_sites)

m = re.search("(?<=n_runs=)\w+", parameters)
n_runs = m.group(0)
n_runs = int(n_runs)

m = re.search("(?<=alphas=).+", parameters)
alphas = m.group(0)
alphas = alphas.split(",")

for j in range(0,len(alphas)):
    alphas[j] = float(alphas[j])

m = re.search("(?<=betas=).+", parameters)
betas = m.group(0)
betas = betas.split(",")

for j in range(0,len(betas)):
    betas[j] = float(betas[j])

m = re.search("(?<=kappa=)\w+", parameters)
kappa = m.group(0)
kappa = int(kappa)

m = re.search("(?<=rate_probs=).+", parameters)
rates = m.group(0)
rates = rates.split(",")

for j in range(0,len(rates)):
    rates[j] = float(rates[j])

############### Loop ##########

for species in species_numbers:
    print species
    check_dir(path.join(os.getcwd(),species))
    os.chdir(path.join(os.getcwd(),species))    
    
    for size in sizes:
        print size
        check_dir(path.join(os.getcwd(),size))
        os.chdir(path.join(os.getcwd(),size))
        
        tree = path.join(treedir, species, size, "tree_file")
        current_tree = pyvolve.read_tree(file = tree)        
        
        for i in range(1,n_runs+1):
            
            os.chdir(path.join(outdir,species,size))
            check_dir(path.join(os.getcwd(),str(i)))
            os.chdir(path.join(os.getcwd(),str(i)))
            
            my_model = pyvolve.Model("codon", {"alpha":alphas, "beta":betas, "kappa":kappa}, rate_probs=rates)
    
            my_partition = pyvolve.Partition(models = my_model, size = n_sites)
            my_evolver = pyvolve.Evolver(partitions = my_partition, tree = current_tree)
            my_evolver()
            







