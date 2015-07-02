from subprocess import Popen
from glob import glob
import os
from os import path
import argparse
import re

import utils

segfault_RE = re.compile("core dumped") # For rerunning

argparser = argparse.ArgumentParser()
argparser.add_argument('--clade', metavar='clade_name', type=str, required=True)
argparser.add_argument('--inroot', metavar='input_root', type=str, required=True)
argparser.add_argument('--treeroot', metavar='tree_root', type=str, required=False, default=None)
argparser.add_argument('--outroot', metavar='out_root', type=str, required=True)
argparser.add_argument('--logdir', metavar='log_dir', type=str, required=True)
argparser.add_argument('--mode', metavar='mode', type=str, required=False, default='codon')
argparser.add_argument('--rerun', action='store_true')

args = argparser.parse_args()

if args.treeroot:
    prank_cmd = "prank -d={} -t={} -o={} -prunetree -" + args.mode
else:
    prank_cmd = "prank -d={} -o={} -prunetree -codon" + args.mode

inroot = args.inroot
treeroot = args.treeroot
alndir = args.outroot
logroot = args.logdir

utils.check_dir(logroot)
utils.check_dir(path.join(logroot, args.clade))
utils.check_dir(alndir)
utils.check_dir(path.join(alndir, args.clade))

for infile in glob(path.join(inroot, args.clade, "*", "*.fa")):
    print infile
    basename = path.basename(infile).partition('.')[0]
    prefix = basename.partition('_')[0][:2]

    outdir = path.join(alndir, args.clade, prefix)
    utils.check_dir(outdir)
    outfile = path.join(outdir, basename + '_prank')

    logdir = path.join(logroot, args.clade, prefix)

    utils.check_dir(logdir)

    logfile = path.join(logdir, basename + '.log')
    errfile = path.join(logdir, args.clade, basename + '.err')

    if treeroot:
        treedir = path.join(treeroot, args.clade, prefix)
        if args.clade == "yeast":
            treefile = path.join(treedir, 'RAxML_bestTree.'+basename)            
        else:
            treefile = path.join(treedir, basename + '.nh')

        prank = prank_cmd.format(infile, treefile, outfile)
    else:
        prank = prank_cmd.format(infile, outfile)

    if args.rerun:
        if segfault_RE.search(open(logfile).read()) is None:
            continue
        else:
            os.remove(logfile)

    p = Popen(['bsub', '-R', 'rusage[tmp=512]', '-o'+logfile, '-g', '/prank',
               '-cwd /tmp', prank])
        
    p.wait()
