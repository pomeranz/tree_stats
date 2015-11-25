from subprocess import Popen
from glob import glob
import os
from os import path
import argparse
import re

import utils
from utils import check_dir

fubar_cmd = "~sparks/hyphy-hyphyqt/HYPHYMP {0}"

argparser = argparse.ArgumentParser()
argparser.add_argument('--inroot', metavar='input_root', type=str, required=True)
argparser.add_argument('--logdir', metavar='log_dir', type=str, required=True)

args = argparser.parse_args()

inroot = args.inroot
logroot = args.logdir

utils.check_dir(logroot)

sizes = "Small", "Medium", "Big"
species_numbers = "6species", "12species", "17species", "44species"

# prepare each of the 12 directories with sequences for slr
for species in species_numbers:
    print species
    check_dir(path.join(inroot,species))    
    
    for size in sizes:
        print size
        check_dir(path.join(inroot, species, size))

        for infile in glob(path.join(inroot, species, size,  "*", "*_FUBAR.ctl")):
            print infile
            
            basename = path.basename(infile)
            logdir = path.join(logroot, basename[:2])
        
            utils.check_dir(logdir)
        
            logfile = path.join(logdir, basename + '.log')
        
            fubar = fubar_cmd.format(infile)
        
            p = Popen(['bsub', '-R', 'rusage[tmp=512]', '-o'+logfile, '-g', '/fubar', fubar])
                
            p.wait()
    
