# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 12:58:30 2015

@author: gideon
"""

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
argparser.add_argument('--outroot', metavar='out_root', type=str, required=True)
argparser.add_argument('--logdir', metavar='log_dir', type=str, required=True)
argparser.add_argument('--rerun', action='store_true')

args = argparser.parse_args()

raxml_cmd = "/homes/pomeranz/standard-RAxML-master/raxmlHPC-PTHREADS-SSE3 -f a -s alg -x 12345 -# 100 -m PROTGAMMAWAGX -s {0} -n {1} -w {2}"

# -n is the name, has to be unique
# -w is the working directory

inroot = args.inroot
outroot = args.outroot
logroot = args.logdir

utils.check_dir(logroot)
utils.check_dir(path.join(logroot, args.clade))
utils.check_dir(outroot)
utils.check_dir(path.join(outroot, args.clade))

for infile in glob(path.join(inroot, "*", "*.fa")):
    print infile
    basename = path.basename(infile).partition('.')[0]
    prefix = basename.partition('_')[0][:2]

    outdir = path.join(outroot, args.clade, prefix)
    utils.check_dir(outdir)
    # outfile = path.join(outdir, basename)

    logdir = path.join(logroot, args.clade, prefix)

    utils.check_dir(logdir)

    logfile = path.join(logdir, basename + '.log')
    errfile = path.join(logdir, args.clade, basename + '.err')

    raxml = raxml_cmd.format(infile, basename, outdir)

    if args.rerun:
        if segfault_RE.search(open(logfile).read()) is None:
            continue
        else:
            os.remove(logfile)

    p = Popen(['bsub', '-R', 'rusage[tmp=512]', '-o'+logfile, '-g', '/raxml',
               '-cwd /tmp', raxml])
        
    p.wait()