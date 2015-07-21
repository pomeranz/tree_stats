# Analyze tdg09 results manually

library(yaml)

# read in the tdg09 output file
out <- yaml.load_file(input='/home/gideon/Documents/mphil_internship/tdg09_test.txt')

# summary 
summary(out)

lrt_results <- as.data.frame(matrix(unlist(out$LrtResults), ncol=5, byrow=T))
names(lrt_results) <- c("site", "deltaLnL", "dof", "lrt", "fdr")

attach(lrt_results)

# fix deltaLnL
deltaLnL = apply(lrt_results, 1, FUN=function(x) abs(x[2]))
lrt_results[2] = deltaLnL

# calculate correct lrt
lrt = apply(lrt_results, 1, FUN=function(x) pchisq(2*x[2], df=x[3], lower.tail=FALSE))
# assign new results
lrt_results[4] = lrt

# calculate4 correct fdr
fdr_vals = apply(lrt_results, 1, FUN=function(x) p.adjust(x[4], method="fdr"))
# assign
lrt_results[5] = fdr_vals

head(lrt_results)

# How many sites did we find with false discovery rate < 0.05?
sum(lrt_results$fdr <= 0.05)

# Prepare the full results table
full_results <- as.data.frame(matrix(unlist(out$FullResults), ncol=9, byrow=T), stringsAsFactors = FALSE)
names(full_results) <- c("site", "ssfParams", "ssfLnL", "lssfParams", "lssfLnL", "deltaLnL", "dof", "lrt", "fdr")

# remove rows with NA
full_results <- full_results[apply(full_results, 1, function(x) {!("NA" %in% x) } ), ]

# make them all to integers
full_results <- apply(full_results, 2, function(x) as.numeric(x))

full_results <- as.data.frame(full_results)

# recalculate same things for this table too
deltaLnL = apply(full_results, 1, FUN=function(x) abs(x[6]))
full_results$deltaLnL = deltaLnL
lrt = apply(full_results, 1, FUN=function(x) pchisq(2*x[6], df=x[7], lower.tail=FALSE))
full_results$lrt = lrt
fdr_vals = apply(full_results, 1, FUN=function(x) p.adjust(x[8], method="fdr"))
full_results$fdr = fdr_vals

# Here's what we have
head(full_results)

# Split the sequence into 5 chunks
sites <- out$Alignment$SiteCount
plot_ranges <- split(seq(1, sites), cut(seq(1, sites), 5))

# Draw the plot
par(mfrow=c(5,1), mar=c(2.0,0.5,0.5,0.5))
for (p in 1:5) {
  plot(1 - fdr,  xlim=c(plot_ranges[[p]][1], tail(plot_ranges[[p]], n=1)), ty='h', lwd=1, main='', xlab='', ylab='', yaxt='n', col="#1B9E77")
  lines(which(fdr <= 0.20), 1 - fdr[fdr <= 0.20], col="#D95F02", ty='h')
  abline(h=0.95, lty='dashed')
  points(which(fdr <= 0.05), 1 - fdr[fdr <= 0.05], pch=20, col="#DE2D26")
}