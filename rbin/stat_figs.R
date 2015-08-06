### Figures for tree_stats.py ###

## PACKAGES ##
library(ggplot2)

#-----------------------------------------------------#

# set wd
setwd("~/Documents/mphil_internship/trees")

# read in the stats file
tree_stats <- read.csv(file = "tree_stats_output.csv", header = TRUE)

colnames(tree_stats)

# make a boxplot of the stats and maybe some histograms of important stuffs

to_remove <- c(1, 5, 7)

tree_numbs <- tree_stats[-to_remove]

attach(tree_numbs)

no_paralogs <- tree_numbs[8]
no_paralogs <- unlist(no_paralogs[1])


# most of the plots are very skewed, need to think of a better way of displaying them. 

# log?
#multiple axis?
# SCALE BREAK, good when we only have a few big numbers. I think that would be best when we display all data at once + multiple axis


## PLOTS ##

## Tree Length ##

plot(Length,type="l", main="Tree Lengths", xlab="Tree #")

## Biggest Branch ##

# maybe x-axis break
hist(Biggest.Branch)

## Smallest Branch ##
hist(Smallest.Branch)

## Distance of farthest leaf ##
plot(Distance.of.farthest.leaf)

## Distance of closest leaf ##
plot(Distance.of.closest.leaf)

## Number of leaves ##
boxplot(X..of.Leaves)

## Number of Species ##
boxplot(X..of.Species)

## Number of Paralogs ##
boxplot(X..of.Paralogs)
hist(log(no_paralogs))

## Number of human seqs ##
boxplot(X..of.Human.seqs)

