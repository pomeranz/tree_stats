import sys
import re 
import os
from os import path
from glob import glob
import pickle
from argparse import ArgumentParser

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import utils

argparser = ArgumentParser()

argparser.add_argument('--clade', metavar='clade_name', type=str, required=True)
argparser.add_argument('--cds', metavar='cds_dir', type=str, required=True)
argparser.add_argument('--pep', metavar='pep_dir', type=str, required=True)
argparser.add_argument('--inroot', metavar='inroot', type=str, required=True)
argparser.add_argument('--outroot', metavar='outroot', type=str, required=True)
argparser.add_argument('--species_cache', metavar='pickle_file', type=str, required=True)
# argparser.add_argument('--specieslist', metavar='species_file', type=str, required=True)

args = argparser.parse_args()

transcript_RE = re.compile("transcript:(\w+)")

# all_species = []
# species_file = args.specieslist
# for l in open(species_file):
#     f = l.strip().split('\t')
#     all_species.append(f[1].lower())
#     print all_species[-1]


print >>sys.stderr, "Building the ID dictionary..."
t2p = {} # Transcript to protein ID
ens_pep_dir = args.pep
for f in glob(path.join(ens_pep_dir, '*.all.fa')):
    print "Processing", f
    for seqr in SeqIO.parse(f, 'fasta'):
        tr = transcript_RE.search(seqr.description)
        if not tr:
            print seqr.description
            raise ValueError("Transcript ID not found in FASTA description field.")

        tr_name = tr.groups()[0]
        t2p[tr_name] = seqr.id

## Build a mapping of Ens IDs to coding sequences 
ens_map = {}
ens_cdna_dir = args.cds
for f in glob(path.join(ens_cdna_dir, '*.fa')):
    print "Processing", f
    for seqr in SeqIO.parse(f, 'fasta'):
        if seqr.id in ens_map:
            print "Duplicate id", seqr.id
            sys.exit(-1)

        ens_map[t2p[seqr.id]] = seqr.seq

clades_pickle = args.species_cache
clades = pickle.load(open(clades_pickle))

inroot = args.inroot
outroot = args.outroot
utils.check_dir(outroot)

for seqset in glob(path.join(inroot, args.clade, "*", "*.tab")):
    setid = path.basename(seqset).rpartition('.')[0]
    seqs = []
    utils.check_dir(path.join(outroot, args.clade))

    for l in open(seqset):
        seqid, species = l.rstrip().split('\t')
        # if species not in all_species:
        #     continue

        # TODO The completely honest thing to do would be to check if the genes
        # with missing sequences fall into the relevant clade or not
        seq = ens_map.get(seqid)
        if seq is None:
            print seqid, "missing!"
            # break
        else:
            seqs.append(SeqRecord(seq, id=seqid, description=""))
    else:
        outdir =  path.join(outroot, args.clade, setid.partition('_')[0][:2])
        utils.check_dir(outdir)

        SeqIO.write(seqs, open(path.join(outdir, setid + '.fa'), 'w'), 'fasta')
