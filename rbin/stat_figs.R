### Figures for tree_stats.py ###

## PACKAGES ##
library(ggplot2)
require(grid)

#-----------------------------------------------------#

# set wd
setwd("~/Documents/mphil_internship/split_tree_out/unprefixed")

# read in the stats file
tree_stats <- read.csv(file = "split_tree_stats.csv", header = TRUE)

colnames(tree_stats)

# make a boxplot of the stats and maybe some histograms of important stuffs

to_remove <- c(1, 5, 7)

tree_numbs <- tree_stats[-to_remove]

attach(tree_numbs)
# Histogram for Length
length <- qplot(Length, xlim=c(0,300))
length + geom_histogram(aes(fill = ..count..))

tail.length <- qplot(X..of.Leaves, xlim=c(300,500), xlab="", ylab = "")
tail.length <- tail.length + geom_histogram(aes(fill = ..count..))
print(tail.length, vp=viewport(width = 0.3, height = 0.3, x = 0.85, y = 0.2, just = c("right", "bottom")))

# Histogram for Number of Leaves
no.leaves <- qplot(X..of.Leaves, data=tree_numbs, xlim=c(0,300), xlab = "Number of leaves")
no.leaves + geom_histogram(aes(fill= ..count..), binwidth=10)

tail.leaves <- qplot(X..of.Leaves, xlim=c(300,500), xlab="", ylab = "")
tail.leaves <- tail.leaves + geom_histogram(aes(fill = ..count..))

print(tail.leaves, vp=viewport(width = 0.3, height = 0.3, x = 0.85, y = 0.2, just = c("right", "bottom")))

# Decided to take this as an example when stuff goes wrong. Have to remove trees where less than 48 leaves. 

# table that holds the means of all the stats

stats_summary <- matrix(ncol=9) 
colnames(stats_summary) <-  c("Length", "Biggest.Branch", "Smallest.Branch", "Farthest.Leaf", "Closest.Leaf", "No.of.Leaves", "No.of.Species", "No.of.Paralogs","No.of.Human.seqs")

for (i in 1:length(colnames(tree_numbs))) {
  stats_summary[1,i] <- round(mean(tree_numbs[,i]))
}
