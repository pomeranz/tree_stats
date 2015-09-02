# Packages
library(yaml)
library(iterators)
library(stringi)
library(stringr)
library(base)
library(xtable)

########################################## TDG09 results extraction function ##################################################
tdg09_results <- function(infile) {
  
  if (file.exists(infile) == FALSE) {
    run = FALSE
    return(list("No file available",run,0))
  }
  
  if (length(readLines(infile)) < 50) {
    run = FALSE
    return(list("No data available",run,0))
  } 
  
  out <- yaml.load_file(input=infile)
  
  if (length(out$FullResults) == 0) {
    run = FALSE
    return(list("Unfinished job",run,0))
  }
  
  # Prepare the full results table
  full_results <- as.data.frame(matrix(unlist(out$FullResults), ncol=9, byrow=T), stringsAsFactors = FALSE)
  names(full_results) <- c("site", "ssfParams", "ssfLnL", "lssfParams", "lssfLnL", "deltaLnL", "dof", "lrt", "fdr")
  
  # make them all to numerics
  full_results <- apply(full_results, 2, FUN=function(x) as.double(x))
  # back to data frame
  full_results <- as.data.frame(full_results)
  
  # recalculate same things for this table too
  deltaLnL_vals = apply(full_results, 1, FUN=function(x) abs(x[6]))
  full_results$deltaLnL = deltaLnL_vals
  #lrt_vals = apply(full_results, 1, FUN=function(x) pchisq(2*x[6], df=x[7], lower.tail=FALSE))
  lrt_vals = apply(full_results, 1, FUN=function(x) pchisq(2*x[6], df=19, lower.tail=FALSE))
  full_results$lrt = lrt_vals
  fdr_vals = p.adjust(full_results[,8], method="fdr")
  full_results$fdr = fdr_vals
  
  fdr_vals[is.na(fdr_vals)] <- 1.0 # conserved locations implicitly have no evidence of non-homogeneity
  
  # Output prep
  signif_vals <- which(fdr_vals <= 0.05)  # Which hits
  no_signifs <- length(signif_vals)  # How many
  
  # Table that has the significant hits
  output <- full_results[signif_vals,]
  
  final_out <- output[,-c(2,3,4,6,7,8)]
  
  run = TRUE
  
  return(list(final_out,run,nrow(full_results)))
  
}

########################################## SLR results extraction function ##################################################

slr_results <- function(infile) {
  
  if (file.exists(infile) == FALSE) {
    run = FALSE
    return(list("No data available", run))
  }
  
  cn <- c("Neutral", "Optimal", "omega", "lower", "upper", "LRT_Stat", "Pval", "Adj.Pval", "Q-value", "Result", "Note")
  
  slr_results <- read.fwf(infile, widths=c(9, 8, 9, 9, 9, 7, 9, 10, 11, 11, 7, 12), header=F, row.names=1, as.is=T)
  
  colnames(slr_results) <- cn
  
  
  signif_vals <- which(slr_results$Adj.Pval <= 0.001 & slr_results$omega > 1.0)  # Which hits
  no_signifs <- length(signif_vals)
  
  output <- slr_results[signif_vals,]
  
  final_out <- output[,-c(1,2,4,5)]
  rownames(final_out) <- gsub(" ", "", substr(rownames(final_out), 1, nchar(rownames(final_out))-1), fixed=TRUE)
  
  run = TRUE
  
  return(list(final_out,run,length(rownames(slr_results))))
}

######################################### RESULTS ##########################################

results_table <- data.frame(matrix(ncol=7))
colnames(results_table) <- c("Alignfile", "TDG09_sites", "PAML_sites", "SLR_sites", "TDG09-PAML", "TDG09-SLR", "SLR-PAML")

setwd("/home/gideon/Documents/mphil_internship/")
inroot = getwd()

# to put results into the right rows. 
index = 1

# Fischers exact test
total_sites <- c()
total_tdg <- c()
total_slr <- c()

# get exact sites of the intersect
exact_intersect <- list()

for (infile in Sys.glob(file.path(inroot, "prank_out", "*/*_prank.best.fas"))) {
  
  print(infile)
  
  # create basenames to acces the other files.
  basename = head(unlist(strsplit(tail(unlist(strsplit(infile, split="/")), n=1), split=".", fixed=TRUE)), n=1)
  basename = paste(unlist(strsplit(basename, split="_"))[1], unlist(strsplit(basename, split="_"))[2], sep="_")
  prefix = substr(unlist(strsplit(basename, split="_", fixed=TRUE))[1], 1, 2)
  
  results_table[index,1] <- paste(basename, ".fa", sep="")
  
  # TDG09
  tdg_file <- paste("/home/gideon/Documents/mphil_internship/tdg09_out/", prefix, "/", basename, ".txt", sep="")
  
  tdg_res_plus <- tdg09_results(tdg_file)
  
  tdg_res <- tdg_res_plus[[1]]
  tdg_sites <- nrow(tdg_res)
  
  if (tdg_res_plus[[2]] == FALSE) {
    results_table[index,2] <- "Missing output"
    tdg_site_names <- c()
    total_tdg <- c(total_tdg, 0)
  } else {
    if (tdg_sites == 0) {
      results_table[index,2] <- tdg_sites
      tdg_site_names <- c()
      total_tdg <- c(total_tdg, 0)
    } else {
      tdg_sites <- nrow(tdg_res)
      results_table[index,2] <- tdg_sites
      tdg_site_names <- tdg_res[,1]
      total_tdg <- c(total_tdg, tdg_sites)
    }
  }

  
  # SLR
  slr_file <- paste("/home/gideon/Documents/mphil_internship/slr_out/Eutheria/",prefix, "/", basename, "_matched.res", sep="")
  
  slr_res_plus <- slr_results(slr_file)
  if (slr_res_plus[[2]] == FALSE) {
    results_table[index,4] <- "Missing output"
    slr_sites <- 0
    slr_site_names <- c()
    total_slr <- c(total_slr, slr_sites)
  } else {
    slr_res <- as.data.frame(slr_res_plus[1])
    
    if (nrow(slr_res) >= 1) {
      slr_sites <- nrow(slr_res)
      results_table[index,4] <- slr_sites
      slr_site_names <- c(as.integer(rownames(slr_res)))
      total_slr <- c(total_slr, slr_sites)
    } else {
      slr_sites <- nrow(slr_res)
      results_table[index,4] <- slr_sites
      slr_site_names <- c()
      total_slr <- c(total_slr, slr_sites)
    }
    
  }
  
  # count number of sites in total and number of sites found in each method
  if (slr_res_plus[[2]] == FALSE) {
    total_sites <- c(total_sites, tdg_res_plus[[3]])
  } else if (length(tdg_sites) == 0) {
    total_sites <- c(total_sites, slr_res_plus[[3]])
  } else if (slr_res_plus[[2]] == FALSE & tdg_res_plus[[2]] == FALSE) {
    total_sites <- c(total_sites, 0)
  }
  
  
  
  # Intersects
  if (length(tdg_site_names) >= 1 & length(slr_site_names) >= 1) {
    get_exact <- TRUE
    tdg_slr_sites <- intersect(tdg_site_names, slr_site_names)
    tdg_slr <- length(intersect(tdg_site_names, slr_site_names))
    if (length(tdg_slr) == 0) {
      results_table[index,6] <- 0
    } else {
      results_table[index,6] <- tdg_slr
      #results_table[index,7] <- tdg_slr_sites
    }
  } else {
    get_exact <- FALSE
    results_table[index,6] <- "Uneven"
  }
  
  # get a list with the exact position of the intersected sites
  if (get_exact == TRUE) {
    title <- basename
    exact_intersect <- c(exact_intersect, list(tdg_slr_sites))
  }
  
  
  get_exact <- FALSE
  
  # progress index
  index = index + 1
  
}

  ######################################### Fisher test ##########################################
# how many times overlapped
total_intersect <- length(which(results_table[,6] != "Uneven" & results_table[,6] != 0))

  
  
m <- matrix(c(total_intersect, sum(total_slr)-total_intersect, sum(total_tdg)-total_intersect, sum(total_sites)-(sum(total_slr)+sum(total_tdg))+total_intersect), nrow=2)
fisher.test(m, alternative="greater")  
  
  
  
