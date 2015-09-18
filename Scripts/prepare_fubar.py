"""
'/home/gideon/Documents/mphil_internship/slr_prep_out/Eutheria/1/1_1_slr.paml'
"""

import os 
from os import path as path
import shutil
import glob
import re
import utils
import argparse
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from dendropy import Tree


# Generate FASTA files
# Copy trees over?
# Make sure file extensions don't clash?

fubar_str = """inputRedirect = {};
inputRedirect["01"]="Universal";
inputRedirect["02"]="1";
inputRedirect["03"]="%s";
inputRedirect["04"]="%s";
inputRedirect["05"]="20";
inputRedirect["06"]="5";
inputRedirect["07"]="2000000";
inputRedirect["08"]="1000000";
inputRedirect["09"]="100";
inputRedirect["10"]="0.5";
ExecuteAFile ("/nfs/research2/goldman/gregs/HBL/FUBAR/FUBAR.bf", inputRedirect);
"""    
argparser = argparse.ArgumentParser()
argparser.add_argument('--indir', metavar='input_directory', type=str, required=True)
argparser.add_argument('--outdir', metavar='input_directory', type=str, required=True)
argparser.add_argument('--clade', metavar='input_directory', type=str, required=True)
args = argparser.parse_args()

utils.check_dir(args.outdir)
utils.check_dir(path.join(args.outdir, args.clade))

def read_slr(fh):
    stats = fh.readline()
    seqs = []

    for l in utils.grouper(fh, 2):
        name = l[0].rstrip()
        seq = l[1].rstrip()
        seqs.append(SeqRecord(id=name, seq=Seq(seq), description=""))
        
    return seqs

for f in glob.glob(path.join(args.indir, args.clade,
                             '*', '*_slr.paml')):
                                 
    dirname, basename = path.split(f)
    input_name = basename.rpartition('.')[0]
    input_core = input_name.rpartition('_')[0]
    seqs = read_slr(open(f))
    
    utils.check_dir(path.join(args.outdir, args.clade, basename[:2]))
    out_fasta = path.abspath(path.join(args.outdir, args.clade, basename[:2], input_core+'.fa'))
    SeqIO.write(seqs, out_fasta, 'fasta')
    
    with open(path.join(dirname, input_name+'.nwk')) as fl:
        _ = fl.readline()
        tree_string = fl.readline().rstrip()
          
    tree = Tree.get_from_string(tree_string, 'newick')    
    
    #tree = Tree.get_from_path(path.join(dirname, input_name+'.nwk'), 'newick', preserve_underscores=True)
    # rewrite the numbered tree back to a tree with ids. 
    # easy
    fasta_file = open(f, 'r')
    fasta = fasta_file.readlines()
    aln_ids = {}
    for index,line in enumerate(fasta,start=1):
        if index%2 == 0:
            aln_ids[(index)/2] = line.strip("\n")
            #print line
            #print (index)/2
        
    for node in tree:
        if node.is_leaf():
            node.taxon.label = aln_ids[int(node.taxon.label)]
            
    # write out the new tree
    out_tree = open(path.abspath(path.join(args.outdir, args.clade, basename[:2], input_core+'.nwk')), 'w')
    out_tree.write(tree.as_string('newick', suppress_rooting=True))

    # Write out the control file
    control_file = open(path.join(args.outdir, args.clade, basename[:2], input_core + '_FUBAR.ctl'), 'w')
    
    print >>control_file, fubar_str % (out_fasta, out_tree)
