from subprocess import Popen
from glob import glob
import os
from os import path
import argparse
import re

import utils

fubar_cmd = "~sparks/hyphy-hyphyqt/HYPHYMP {0}"

argparser = argparse.ArgumentParser()
argparser.add_argument('--clade', metavar='clade_name', type=str, required=True)
argparser.add_argument('--inroot', metavar='input_root', type=str, required=True)
argparser.add_argument('--logdir', metavar='log_dir', type=str, required=True)

args = argparser.parse_args()

inroot = args.inroot
logroot = args.logdir

utils.check_dir(logroot)
utils.check_dir(path.join(logroot, args.clade))

for infile in glob(path.join(inroot, args.clade, "*", "*_FUBAR.ctl")):
    print infile
    
    basename = path.basename(infile)
    logdir = path.join(logroot, args.clade, basename[:2])

    utils.check_dir(logdir)

    logfile = path.join(logdir, basename + '.log')

    fubar = fubar_cmd.format(infile)

    p = Popen(['bsub', '-R', 'rusage[tmp=512]', '-o'+logfile, '-g', '/fubar', fubar])
        
    p.wait()
    
