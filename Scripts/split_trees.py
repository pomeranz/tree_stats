import os
from os import path
import sys
import pickle
import operator
import time
from pprint import pprint
from argparse import ArgumentParser

import ete2
import emf
import utils

argparser = ArgumentParser()

# EMF file is /nfs/research2/goldman/gregs/slr_pipeline/data/ens/78/Compara.78.protein.nhx
# Clade is Eutheria
# species_clade is
# /nfs/research2/goldman/gregs/slr_pipeline/data/ens/78/clades.pk
argparser.add_argument('--emf', metavar='emf_file', type=str, required=True)
argparser.add_argument('--clade', metavar='clade_name', type=str, required=True)
# argparser.add_argument('--treeroot', metavar='tree_root', type=str, required=True)
# argparser.add_argument('--imgroot', metavar='img_root', type=str, required=True)
argparser.add_argument('--outroot', metavar='out_root', type=str, required=True)
argparser.add_argument('--species_cache', metavar='cache_file', type=str, required=True)
argparser.add_argument('--thr', metavar='value', type=float, default=0.6)
# add argument if you want prefixes 
argparser.add_argument('--prefix', metavar='add prefix?', type=str, required=True)


def filter_clades(species, clade_names):
    clade_names = set(clade_names)
    clades = {}
    for cl in clade_names:
        clades[cl] = []
        
    parent = None
    for sp in species:
        sp_clades = clade_names.intersection(get_clades(sp))
        for cl in sp_clades:
            clades[cl].append(sp)

    return clades

clade_cache = {}
def get_clades(sp):
    parent = sp
    clades = []
    while parent != "33154": # Eukaryota
        if parent not in clade_cache:
            if parent == "canis familiaris":
                parent = "canis lupus familiaris"
            sp_info = utils.ens_get("taxonomy/id/", parent)
            clade_cache[parent] = sp_info
            time.sleep(0.33)
        else:
            print parent, "cached!"
            sp_info = clade_cache[parent]
        # pprint(sp_info["parent"])
        clades.append(sp_info["name"])
        parent = sp_info["parent"]["id"]

    return clades

def load_clade(fh):
    for l in fh:
        name, taxid = l.rstrip().split('\t')
        self.thr = thr

class TCC(object):
    def __init__(self, clade, condition, thr):
        self.clade = clade
        self.condition = condition
        self.thr = thr

class TCCList(object):
    def __init__(self):
        self.tccs = []

    def add(self, tcc):
        self.tccs.append(tcc)

    def check(self, node):
        for tcc in self.tccs:
            desc_species = set([ n.species for n in node.get_leaves() ])
            clade_cov = float(len(desc_species.intersection(tcc.clade))) / len(tcc.clade)
            if not tcc.condition(clade_cov, tcc.thr):
                return False

        return True

def split_tree(intree, tcclist, prefix):
    seqsets, subtrees = [], []

    trees = []
    for node in intree.traverse("postorder"):
        if node.dist > 100:
            node.detach()
            if len(node.get_descendants()) > 0:
                trees.append(node)
        # TODO: Does this do anything?
        elif node.is_leaf() and not getattr(node, 'species', None):
            print "REMOVING b/c of lack of species ann."
            node.detach()

    if len(intree.get_children()):
        trees.append(intree)
    # if len(to_remove):
    #     print to_remove
    #     for n in to_remove:
    #         n.detach()

    # TMP
    for node in intree.get_leaves():
        if not getattr(node, 'species', None):
            print "MISSING SPECIES", node.name
            print node.get_children()
    

    for tree in trees:
        
        for node in tree.traverse("postorder"):
            children = node.get_children()
            # changed any to all. 
            if any([ n.done for n in children ]):
                # if not all([ n.done for n in children ]):
                #     print "Warning!"
                node.add_features(done=True)
                continue
             ##rework for if all the childs meet criteria add node # DONE
            if len(children) and all([ n.split for n in children ]):
                                
                seqsets.append([ (n.name, n.species) for n in node.get_leaves() ])
                subtrees.append(node)

                node.add_features(done=True)
            else:
                if tcclist.check(node):
                    node.add_features(split=True)
                else:
                    node.add_features(split=False)

                node.add_features(done=False)
                
    # add the prefix for the particular tree
                # this method doesnt work...
                # node.get_leaf_names = [prefix + leafname for leafname in node.get_leaf_names()]
                #for identifier in node.traverse():
                    #identifier.name = "".join([prefix + identifier.name])
    
    if prefix == "True":            
        index = 1
        for subtree in subtrees:
            if index == len(subtrees):
                for identifier in subtree.traverse():
                    prefix = "".join(["D" + str(index) + "_"])
                    identifier.name = "".join([prefix + identifier.name])
                index = 1
            else:
                for identifier in subtree.traverse():
                    prefix = "".join(["D" + str(index) + "_"])
                    identifier.name = "".join([prefix + identifier.name])
            index += 1


    return seqsets, subtrees

colours = { -1: "black", 0: "red", 1: "green", 2: "blue", 3: "purple", 4: "orange", 5: "yellow", 6: "grey", 7: "#009999",
            8: "#FF7400"}
def make_layout(nodesets):
    nodemap = {}
    for i, ns in enumerate(nodesets):
        for n in ns:
            nodemap[n[0]] = i

    def layout(node):
        if node.D == "Y":
            node.img_style["shape"] = "circle"
            node.img_style["size"] = 12
            node.img_style["fgcolor"] = "red"

        if node.is_leaf():
            nameFace = ete2.faces.AttrFace("name", fsize=20, 
                                           fgcolor=colours[nodemap.get(node.name, -1)])
            ete2.faces.add_face_to_node(nameFace, node, column=0)

    return layout

def main():
    args = argparser.parse_args()

    emf_file = args.emf
    out_root = args.outroot
    # img_root = args.imgroot
    # tree_root = args.treeroot
    clades_pickle = args.species_cache
    prefix = args.prefix

    all_species = utils.ens_get("/info/species/")["species"]
    all_species_names = [ it["name"].replace("_", " ") for it in all_species ]
    all_species_names.remove("Ancestral sequences")

    if path.exists(clades_pickle):
        Clades = pickle.load(open(clades_pickle, 'rb'))
    else:
        Clades = filter_clades(all_species_names,
                               [ "Eutheria", "Glires", "Laurasiatheria", "Sauria", "Mammalia", "Primates" ])
        pickle.dump(Clades, open(clades_pickle, 'wb'))

    pprint(Clades)

    TL = TCCList()
    TL.add(TCC(Clades[args.clade], operator.ge, args.thr))

    utils.check_dir(path.join(out_root, args.clade))

    tree_id = 1
    for tree in emf.EMF(emf_file):
        print tree_id
        # treedir = path.join(tree_root, str(tree_id)[:2])
        # utils.check_dir(treedir)

        # tree.write(outfile=path.join(treedir, "{}.nh".format(tree_id)))

        seqsets, subtrees = split_tree(tree, TL, prefix)
        outdir = path.join(out_root, args.clade, str(tree_id)[:2])
        utils.check_dir(outdir)

        # Treevis
        # layout = make_layout(seqsets)
        # imgdir = path.join(img_root, args.clade)
        # utils.check_dir(imgdir)
        # imgfile = path.join(imgdir, "{}.pdf".format(tree_id))
        # tree.render(imgfile, layout=layout)

        set_id = 1
        for seqset, subtree in zip(seqsets, subtrees):
            outfile = open(path.join(outdir, "{0}_{1}.tab".format(tree_id, set_id)), 'w')
            for seqid in seqset:
                print >>outfile, '\t'.join(seqid)

            subtree.write(outfile=path.join(outdir, "{0}_{1}.nh".format(tree_id, set_id)), 
                          format=6)
            set_id += 1

        tree_id += 1

if __name__ == "__main__":
    main()
