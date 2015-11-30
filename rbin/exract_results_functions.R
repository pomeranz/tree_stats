### exract_results_functions.R ###

# Script that contains the functions that extract the results from each method. 

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
    return("No file available")
  }
  
  if (length(readLines(infile)) < 50) {
    return("No data available")
  } 
  
  out <- yaml.load_file(input=infile)
  
  if (length(out$FullResults) == 0) {
    return("Unfinished job")
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
  
  return(list(final_out,nrow(full_results)))
  
}


########################################## SLR results extraction function ##################################################

slr_results <- function(infile) {
  
  if (file.exists(infile) == FALSE) {
    return("No data available")
  }
  
  cn <- c("Neutral", "Optimal", "omega", "lower", "upper", "LRT_Stat", "Pval", "Adj.Pval", "Q-value", "Result", "Note")
  
  slr_results <- read.fwf(infile, widths=c(9, 8, 9, 9, 9, 7, 9, 10, 11, 11, 7, 12), header=F, row.names=1, as.is=T)

  colnames(slr_results) <- cn
  
  
  signif_vals <- which(slr_results$Adj.Pval <= 0.001 & slr_results$omega > 1.0)  # Which hits
  no_signifs <- length(signif_vals)
  
  output <- slr_results[signif_vals,]
  
  final_out <- output[,-c(1,2,4,5)]
  rownames(final_out) <- gsub(" ", "", substr(rownames(final_out), 1, nchar(rownames(final_out))-1), fixed=TRUE)
  
  return(list(final_out,length(rownames(slr_results))))
}

########################################## PAML results extraction function ##################################################

paml_results <- function(infile) {
  
  if (file.exists(infile) == FALSE) {
    return("No data available")
  }
  
  if (length(readLines(infile)) < 50) {
    return("No data available")
  }
  
  index = 1
  check = c()
  
  for (i in readLines(infile)) {
    if ("Positively selected sites" == i) {
      check[index] <- "True"
    } else {
      check[index] <- "False"
    }
    index = index + 1
  }
  
  if (!(any(check == "True"))) {
    print("No output")
  }
  
  # paml_res <- read.table(infile, col.names = paste0("V",seq_len(23)), fill=TRUE, skip=1500, header=FALSE, sep=" ")
  
  df <- data.frame(matrix(ncol=3))
  index = 1
  BEB <- FALSE
  go <- FALSE
  
  for (i in stri_read_lines(infile)) {
    
    if ("(BEB)" %in% unlist(str_split(i, pattern= " "))) {
      BEB <- TRUE
    }
    
    if ("Positively" %in% unlist(str_split(i, pattern= " ")) && BEB == TRUE) {
      #print(TRUE)
      #print(unlist(str_split(i, pattern= " ")))
      #temp = i
      go <- TRUE
    }
    
    if (go == TRUE) {
      if (stri_startswith_fixed(i, pattern="  ")) {
        #print(i)
        df[index,] <- unlist(str_split(i, pattern= "     "))
        index = index + 1
      }
    }
  
  # for loop end
  }
  
  colnames(df) <- c("site", "Prob(w>1)", "mean w")
  
  
  
  return(df)
# function end  
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
total_intersect <- list()


for (infile in Sys.glob(file.path(inroot, "alignments/fasta_AA", "*/*.fa"))) {
  
  print(infile)
  
  # create basenames to acces the other files.
  basename = head(unlist(strsplit(tail(unlist(strsplit(infile, split="/")), n=1), split=".", fixed=TRUE)), n=1)
  prefix = substr(unlist(strsplit(basename, split="_", fixed=TRUE))[1], 1, 2)
  
  results_table[index,1] <- paste(basename, ".fa", sep="")
  
  # TDG09
  tdg_file <- paste("/home/gideon/Documents/mphil_internship/tdg09_out/", prefix, "/", basename, ".txt", sep="")
  
  tdg_res_plus <- tdg09_results(tdg_file)
  
  tdg_res <- tdg_res_plus[[1]]
  
  tdg_sites <- nrow(tdg_res)
  if (length(tdg_sites) == 0) {
    no_data <- "Missing output"
    results_table[index,2] <- no_data
  } else if (tdg_sites == 0) {
    results_table[index,2] <- tdg_sites
  } else {
      tdg_sites <- nrow(tdg_res)
      results_table[index,2] <- tdg_sites
  }
  print("TDG done")
  
  # PAML
  paml_file <- paste("/home/gideon/Documents/mphil_internship/paml_out/out/",prefix, "/", basename, "/", "rst", sep="")
  
  paml_res <- paml_results(paml_file)
  
  # this does not return a logical value???
  if (paml_res == "No data available") {
    no_data <- "Missing output"
    results_table[index,3] <- no_data
  } else if (paml_res == "Unfinished job") {
    no_out <- "Unfinished job"
    results_table[index,3] <- no_data
  } else {
    paml_sites <- nrow(paml_res)
    results_table[index,3] <- paml_sites
  }
  print("PAML done")
  
  # SLR
  slr_file <- paste("/home/gideon/Documents/mphil_internship/slr_out/Eutheria/",prefix, "/", basename, "_matched.res", sep="")
  
  slr_res_plus <- slr_results(slr_file)
  slr_res <- as.data.frame(slr_res_plus[1])
  
  if (length(slr_res) > 1) {
    slr_sites <- nrow(slr_res)
    results_table[index,4] <- slr_sites
  } else if (slr_res == "No data available") {
    no_data <- "Missing output"
    results_table[index,4] <- no_data
  } else {
    slr_sites <- nrow(slr_res)
    results_table[index,4] <- slr_sites
  }
  print("SLR done")

  # Calculate fisher tests for TDG-SLR
#   if (tdg_res_plus == "No data available") {
#     N <- slr_res_plus[[2]]
#   } else {
#     N <- tdg_res_plus[[2]] # total number of items
#   }
#   total_sites <- c(total_sites, N)
#   n <- tdg_sites # total white
#   total_tdg <- c(total_tdg, n)
#   k <- slr_sites # total drawn
#   total_slr <- c(total_slr, k)
#   
#   if (tdg_res == "No data available") {
#     TSGs <- 0
#   } else {
#     TDGs <- as.vector(tdg_res$site)
#   } 
#   SLRs <- as.vector(rownames(slr_res))
#   x <- intersect(TDGs, SLRs)
#   total_intersects <- c(total_intersects, x)
#   
#   if (class(x) == "character") {
#     x <- length(x)
#   }
#   results_table[index, 6] <- x
#   
#   print("Fischers done")
  # increase index
  index = index + 1
  
# end of for loop  
}

# calculate fishers exact
# m <- matrix(c(x, k-x, n-x, N-(k+n)+x), nrow=2)
# fisher.test(m, alternative="greater")
# 
# 
# #write.table(results_table, file="/home/gideon/Documents/mphil_internship/results_table.csv", row.names=FALSE)
# 
#sub_results <- xtable(head(results_table[,1:4], n =20))
# results_table[,which(results_table$PAML_sites != "Missing output")]
#sub_results2 <- xtable(results_table[,1:4])

#print(toLatex(sub_results))


