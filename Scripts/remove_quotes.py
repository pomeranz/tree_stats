# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 16:31:04 2015

@author: gideon

Description: This script removes all quotes from the raxml trees.
"""

from glob import glob
import os
from os import path
import argparse


argparser = argparse.ArgumentParser()

argparser.add_argument('--treeroot', metavar='tree_root', type=str, required=True)

args = argparser.parse_args()

for infile in glob(path.join(args.treeroot, "*", "*.nwk")):
    
    tree_file = open(infile, "r+")
    
    text = tree_file.read()
    text = text.replace("'", "")
    tree_file.seek(0)
    tree_file.write(text)
    tree_file.truncate()
    tree_file.close()
