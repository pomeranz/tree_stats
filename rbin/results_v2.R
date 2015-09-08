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
    return(list("No data available", run, 0))
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

########################################## PAML results extraction function ################

# Greg stuff
# grep("(BEB)", readLines(paml_file))
# which(readLines(paml_file) == "")
# beb <- grep("(BEB)", readLines(paml_file))+3
# empty <- which(readLines(paml_file) == "")
# beb_end <- min(empty[empty > beb])
# readLines(paml_file)[beb:beb_end]

paml_results <- function(infile) {
  
  # Error 1 - Didn't run --> no file
  if (file.exists(infile) == FALSE) {
    run = FALSE
    return(list("Error: No outfile", run, 0))
  }
  
  # Error 2 - Didn't start --> less than 50 lines
  if (length(readLines(infile)) < 50) {
    run = FALSE
    return(list("Error: no output", run, 0))
  }
  
  file <- readLines(paml_file)
  
  # Error 3 - unfinished output --> job was killed half way through
  if (!("(BEB)" %in% unlist(str_split(file, pattern= " ")))) {
    run = FALSE
    return(list("Error: unfinished job", run))
  }
  
  
  # locating correct part of the file
  beb <- grep("(BEB)", file)
  end_of_file <- length(file)
  sub_file <- file[beb:end_of_file]
  end_of_subfile <- length(sub_file)
  positive <- grep("Positively", sub_file)+4
  
  # setting up loop
  index = 1
  df <- data.frame(matrix(ncol=3))
  
  for (i in sub_file[positive:end_of_subfile]) {
    if (length(unlist(str_split(i, pattern= " "))) > 30) {
      break
    } else {
      site <- as.numeric(substr(i, 1,7))
      prob <- as.numeric(substr(i, 14,22))
      omega <- substr(i,29,43)
      df[index,] <- c(site,prob,omega)
      index = index + 1
    }
  }
  

  
  colnames(df) <- c("site", "Prob(w>1)", "mean w")
  
  run = TRUE
  
  return(list(df, run, nrow(df)))
  # function end  
}


######################################### RESULTS FUNCTION #########################################################

site_table_rows <- function(file, method, sites, site_table) {
  
  # create basenames to acces the other files.
  basename = head(unlist(strsplit(tail(unlist(strsplit(infile, split="/")), n=1), split=".", fixed=TRUE)), n=1)
  basename = paste(unlist(strsplit(basename, split="_"))[1], unlist(strsplit(basename, split="_"))[2], sep="_")
  prefix = substr(unlist(strsplit(basename, split="_", fixed=TRUE))[1], 1, 2)
  
  if (method == "TDG09") {
    file <- paste("/tdg09_out/", prefix, "/", basename, ".txt", sep="")
  } else if (method == "SLR") {
    file <- paste("/slr_out/Eutheria/",prefix, "/", basename, "_matched.res", sep="")
  } else if (method == "PAML") {
    file <- paste("/paml_out/out/",prefix, "/", basename, "/", "rst", sep="")
  }
  
  if (length(sites) >= 1) {
    for (site in sites) {
      add <- c(file,method,site)
      site_table <- rbind(site_table, add)
    }
  }
  
  return(site_table)
}

######################################### FISHER TEST AND INTERSECT FUNCTIONS ##########################################

count_total_sites <- function(results1,results2, total_sites) {
  
  if (results2[[2]] == FALSE) {
    total_sites <- c(total_sites, results1[[3]])
  } else if (results1[[2]] == FALSE) {
    total_sites <- c(total_sites, results2[[3]])
  } else if (results2[[2]] == FALSE & results1[[2]] == FALSE) {
    total_sites <- c(total_sites, 0)
  }
  
  return(total_sites)
}

# I don't think this is working
site_intersects <- function(method1_sites, method2_sites, results_table, index, column) {
  # calculates the intersect bwtween two methods and adds it to the results_table
  
  # Intersects
  if (length(method1_sites) >= 1 & length(method2_sites) >= 1) {
    get_exact <- TRUE
    intersected_sites <- intersect(method1_sites, method2_sites)
    no_intersected_sites <- length(intersect(method1_sites, method2_sites))
    if (length(no_intersected_sites) == 0) {
      results_table[index,column] <- 0
    } else {
      results_table[index,column] <- no_intersected_sites
      #results_table[index,7] <- tdg_slr_sites
    }
  } else {
    get_exact <- FALSE
    results_table[index,column] <- "Uneven"
  }
  
  return(list(results_table,get_exact))
  
}

methods_fishers <- function(results_table,column,toal_method1,total_method2, total_sites) {
  
  # how many times overlapped
  total_intersect <- length(which(results_table[,column] != "Uneven" & results_table[,column] != 0))
  
  m <- matrix(c(total_intersect, sum(total_method2)-total_intersect, sum(toal_method1)-total_intersect, sum(total_sites)-(sum(total_method2)+sum(toal_method1))+total_intersect), nrow=2)
  x <- fisher.test(m, alternative="greater")
  
  return(x)
}


######################################### RESULTS ##########################################

### Preparation ###
results_table <- data.frame(matrix(ncol=7))
colnames(results_table) <- c("Alignfile", "TDG09_sites", "PAML_sites", "SLR_sites", "TDG09-PAML", "TDG09-SLR", "SLR-PAML")
site_table <- data.frame(matrix(ncol=3))
colnames(site_table) <- c("file", "method", "site")

setwd("/home/gideon/Documents/mphil_internship/")
inroot = getwd()

# to put results into the right rows. 
index = 1

# Fischers exact test
total_sites_tdg_slr <- c()
total_sites_slr_paml <- c()
total_sites_tdg_paml <- c()
total_tdg <- c()
total_slr <- c()
total_paml <- c()

# get exact sites of the intersect
exact_intersect <- list()
intersect_names <- c()

for (infile in Sys.glob(file.path(inroot, "prank_out", "*/*_prank.best.fas"))) {
  
  print(infile)
  
  # create basenames to acces the other files.
  basename = head(unlist(strsplit(tail(unlist(strsplit(infile, split="/")), n=1), split=".", fixed=TRUE)), n=1)
  basename = paste(unlist(strsplit(basename, split="_"))[1], unlist(strsplit(basename, split="_"))[2], sep="_")
  prefix = substr(unlist(strsplit(basename, split="_", fixed=TRUE))[1], 1, 2)
  
  results_table[index,1] <- paste(basename, ".fa", sep="")
  
  ### TDG09 ###
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
  
  # add to site_table
  site_table <- site_table_rows(infile, "TDG09", tdg_site_names, site_table)
  
  ### SLR ###
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
  
  # add to site_table
  site_table <- site_table_rows(infile, "SLR", slr_site_names, site_table)
  
  ### PAML ###
  paml_file <- paste("/home/gideon/Documents/mphil_internship/paml_out/out/",prefix, "/", basename, "/", "rst", sep="")
  
  paml_res_plus <- paml_results(paml_file)
  
  if (paml_res_plus[[2]] == FALSE) {
    results_table[index,3] <- "Missing output"
    paml_sites <- 0
    paml_site_names <- c()
    total_paml <- c(total_paml, paml_sites)
  } else {
    paml_res <- as.data.frame(paml_res_plus[[1]])
    
    if (nrow(paml_res) >= 1) {
      paml_sites <- nrow(paml_res)
      results_table[index,3] <- paml_sites
      paml_site_names <- c(as.integer(rownames(paml_res)))
      total_paml <- c(total_paml,paml_sites)
    } else {
      paml_sites <- nrow(paml_res)
      results_table[index,3] <- paml_sites
      paml_site_names <- c()
      total_paml <- c(total_paml,paml_sites)
    }
  }
  
  # add to site_table
  site_table <- site_table_rows(infile, "PAML", paml_site_names, site_table)
  
  ###### FISHER TESTS AND INTERSECTS #######
  # count number of sites in total and number of sites found in each method
#   if (slr_res_plus[[2]] == FALSE) {
#     total_sites <- c(total_sites, tdg_res_plus[[3]])
#   } else if (length(tdg_sites) == 0) {
#     total_sites <- c(total_sites, slr_res_plus[[3]])
#   } else if (slr_res_plus[[2]] == FALSE & tdg_res_plus[[2]] == FALSE) {
#     total_sites <- c(total_sites, 0)
#   }
  
  # TDG - SLR 
  total_sites_tdg_slr <- count_total_sites(tdg_res_plus, slr_res_plus, total_sites_tdg_slr)
  tdg_slr_intersects <- site_intersects(tdg_site_names, slr_site_names, results_table, index, 6)

  # TDG - PAML
  total_sites_tdg_paml <- count_total_sites(tdg_res_plus, paml_res_plus, total_sites_tdg_paml)
  tdg_paml_intersects <- site_intersects(tdg_site_names, paml_site_names, results_table, index, 5)

  # SLR - PAML
  total_sites_slr_paml <- count_total_sites(slr_res_plus, paml_res_plus, total_sites_slr_paml)
  slr_paml_intersects <- site_intersects(slr_site_names, paml_site_names, results_table, index, 7)
  
  
#   # Intersects
#   if (length(tdg_site_names) >= 1 & length(slr_site_names) >= 1) {
#     get_exact <- TRUE
#     tdg_slr_sites <- intersect(tdg_site_names, slr_site_names)
#     tdg_slr <- length(intersect(tdg_site_names, slr_site_names))
#     if (length(tdg_slr) == 0) {
#       results_table[index,6] <- 0
#     } else {
#       results_table[index,6] <- tdg_slr
#       #results_table[index,7] <- tdg_slr_sites
#     }
#   } else {
#     get_exact <- FALSE
#     results_table[index,6] <- "Uneven"
#   }

  
  # get a list with the exact position of the intersected sites
#   if (get_exact == TRUE) {
#     if (length(tdg_slr_sites) > 1) {
#       intersect_names <- c(intersect_names, basename) # not tested yet
#       exact_intersect <- c(exact_intersect, list(tdg_slr_sites))
#     } else {
#       exact_intersect <- c(exact_intersect, tdg_slr_sites)
#     }
#     
#   }
  
  
  #get_exact <- FALSE
  
  # progress index
  index = index + 1
  
}

######################################### Fisher test ##########################################
# how many times overlapped
# total_intersect <- length(which(results_table[,6] != "Uneven" & results_table[,6] != 0))
#   
# m <- matrix(c(total_intersect, sum(total_slr)-total_intersect, sum(total_tdg)-total_intersect, sum(total_sites)-(sum(total_slr)+sum(total_tdg))+total_intersect), nrow=2)
# fisher.test(m, alternative="greater") 

## TDG-SLR ##
methods_fishers(results_table, 6, total_tdg, total_slr, sum(total_sites_tdg_slr))

## TDG-PAML ##
methods_fishers(results_table, 6, total_tdg, total_paml, sum(total_sites_tdg_paml))

## SLR-PAML ##
methods_fishers(results_table, 6, total_slr, total_paml, sum(total_sites_slr_paml))
  
######################################### Intersect shenanigans ########################################## 
names(exact_intersect) <- intersect_names
exact_intersect

######################################### Inspect results ########################################## 

# results_table
site_table <- site_table[-1,]
#site_table
