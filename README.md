# tree_stats
A script that takes a number of tree directories and outputs tree statistics

When the script is run it creates a function tree_stats. 

It needs to be parsed the argument "directory", which tells it what directory the tree files are in.
The two optional arguments are "outname" and "outdir", which specify the filename and directroy of the output file, respectively. 

The tabular file can than be found in the specified directory, or if not in the current working directory. 

While running, the script prints the current tree working directory.

How to run:

python tree_stats --directory "Your tree directory" [--outname] [--outdir]
