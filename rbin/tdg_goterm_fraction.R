# Packages
library(yaml)
library(biomaRt)
library(iterators)
library(stringi)
library(stringr)
library(base)
library(ape)
library(seqinr)

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

#############################################################################################################

### Preparation ###
results_table <- data.frame(matrix(ncol=5))
colnames(results_table) <- c("Alignfile", "Alignment_Length", "TDG09_sites","len_isect", "goTerm_frac")
site_table <- data.frame(matrix(ncol=3))
colnames(site_table) <- c("file", "method", "site")

# set directory right
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

# set up biomart
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="www.ensembl.org", path="/biomart/martservice", dataset = "hsapiens_gene_ensembl")
evidence_types <- c("EXP", "IDA", "IPI", "IMP", "IGI", "IEP")

# setup for go term comparison
found_sites <- c()
no_sites <- c()

# write a loop where it simply records all the human ids in the entire dataset

all_human_ids <- data.frame(matrix(ncol=3))
colnames(all_human_ids) <- c("infile", "D_id", "id")

index = 1

for (infile in Sys.glob(file.path(inroot, "alignments", "fasta_AA", "*/*.fa"))) {
  
  print(infile)
  
  # create basenames to acces the other files.
  basename = head(unlist(strsplit(tail(unlist(strsplit(infile, split="/")), n=1), split=".", fixed=TRUE)), n=1)
  basename = paste(unlist(strsplit(basename, split="_"))[1], unlist(strsplit(basename, split="_"))[2], sep="_")
  prefix = substr(unlist(strsplit(basename, split="_", fixed=TRUE))[1], 1, 2)
  
  # get human ids
  fasta <- read.fasta(file=infile)
  seq_ids <- names(fasta)
  human_ids <- seq_ids[grepl("ENSP0", seq_ids)]
  
  
  for (i in 1:length(human_ids)) {
    D_id <- human_ids[i]
    id <- unlist(strsplit(D_id, "_"))[2]
    all_human_ids <- rbind(all_human_ids, c(infile,D_id,id))
  }
  
  
  index = index + 1
}

# remove first row
all_human_ids <- all_human_ids[2:nrow(all_human_ids),]
rownames(all_human_ids) <- as.integer(rownames(all_human_ids)) - 1

# try bimart?
goterms <- getBM(attributes = c("ensembl_peptide_id", "go_id", "name_1006", "definition_1006", "go_linkage_type", "namespace_1003"),
                    filters = list(ensembl_peptide_id=all_human_ids$id),
                    mart = mart)

# reset the index
index = 1

# now loop through the alignments again and access the goterm data from the goterms dataframe. 
for (infile in Sys.glob(file.path(inroot, "alignments", "fasta_AA", "*/*.fa"))) {
  
  print(infile)
  
  # get alignment length
  aln_length <- as.integer(summary(read.FASTA(infile))[,1][1])
  results_table[index,2] <- aln_length
  
  # create basenames to acces the other files.
  basename = head(unlist(strsplit(tail(unlist(strsplit(infile, split="/")), n=1), split=".", fixed=TRUE)), n=1)
  basename = paste(unlist(strsplit(basename, split="_"))[1], unlist(strsplit(basename, split="_"))[2], sep="_")
  prefix = substr(unlist(strsplit(basename, split="_", fixed=TRUE))[1], 1, 2)
  
  results_table[index,1] <- paste(basename, ".fa", sep="")
  
  ### TDG09 ###
  tdg_file <- paste("/home/gideon/Documents/mphil_internship/tdg09/tdg09_out/", prefix, "/", basename, ".txt", sep="")
  
  tdg_res_plus <- tdg09_results(tdg_file)
  
  tdg_res <- tdg_res_plus[[1]]
  tdg_sites <- nrow(tdg_res)
  
  
  # GO terms business
  
  fasta <- read.fasta(file=infile)
  seq_ids <- names(fasta)
  human_ids <- seq_ids[grepl("ENSP0", seq_ids)]
  
  # get D1 human ids
  D1 <- human_ids[grepl("D1_", human_ids)]
  D1_ids <- c()
  
  for (i in strsplit(D1, "_")) {
    D1_ids <- c(D1_ids, i[2])
  }
  
  if (length(D1_ids) == 0) {
    fract = 0
    next
  } else {
    terms_D1 <- c()
    
    for (j in D1_ids) {
      term_D1 <- goterms[which(goterms$ensembl_peptide_id == j),]
      term_D1 <- subset(term_D1, namespace_1003 == "molecular_function")
      
      terms_D1 <- rbind(terms_D1, term_D1)
    }
    
    #sum(results_D1$go_linkage_type %in% evidence_types)
    
    results_subset_D1 <- subset(terms_D1, go_linkage_type %in% evidence_types)
  }
  
  
  # get D2 human ids
  
  D2 <- human_ids[grepl("D2_", human_ids)]
  D2_ids <- c()
  
  for (i in strsplit(D2, "_")) {
    D2_ids <- c(D2_ids, i[2])
  }
  
  if (length(D2_ids) == 0) {
    fract = 0
    next
  } else {
    
    terms_D2 <- c()
    
    for (j in D2_ids) {
      term_D2 <- goterms[which(goterms$ensembl_peptide_id == j),]
      term_D2 <- subset(term_D2, namespace_1003 == "molecular_function")
      
      terms_D2 <- rbind(terms_D2, term_D2)
    }
    
    #sum(results_D2$go_linkage_type %in% evidence_types)
    
    results_subset_D2 <- subset(terms_D2, go_linkage_type %in% evidence_types)
  }
  
  if (length(results_subset_D1$name_1006) == 0 | length(results_subset_D2$name_1006) == 0) {
    isect=NA
    fract=NA
  } else {
    # intersect stuuf
    isect <- intersect(results_subset_D1$name_1006,results_subset_D2$name_1006)
    # calculate the fraction
    fract <- as.numeric(length(isect))/length(union(results_subset_D1$name_1006,results_subset_D2$name_1006))
  }
    
  if (length(fract) == 0) {
    fract = NA
  }
  
  
  
  if (tdg_res_plus[[2]] == FALSE) {
    results_table[index,3] <- NA
    tdg_site_names <- c()
    total_tdg <- c(total_tdg, 0)
    results_table[index,4] <- NA
    results_table[index,5] <- NA
  } else {
    if (tdg_sites == 0) {
      results_table[index,3] <- tdg_sites
      tdg_site_names <- c()
      total_tdg <- c(total_tdg, 0)
      no_sites <- c(no_sites, fract)
      results_table[index,4] <- length(isect)
      results_table[index,5] <- fract
    } else {
      tdg_sites <- nrow(tdg_res)
      results_table[index,3] <- tdg_sites
      tdg_site_names <- tdg_res[,1]
      total_tdg <- c(total_tdg, tdg_sites)
      found_sites <- c(found_sites,fract)
      results_table[index,4] <- length(isect)
      results_table[index,5] <- fract
    }
  }
  
  # add to site_table
  site_table <- site_table_rows(infile, "TDG09", tdg_site_names, site_table)
  
  index <- index + 1
}

################################ ANALYSIS ################################

# for now no_sites is holding the ones where no_sites where found and 
# found_sites shows the ones that did. 

# count the number of times each number occurs. 

plot(density(table(no_sites)), col="blue", lty=2)
lines(density(table(found_sites)), col="red")

# remove any NA

is.na(results_table)
results_subset <- subset(results_table, !is.na(results_table$goTerm_frac))
with_sites <- subset(results_subset, results_subset$TDG09_sites >= 1)
without_sites <- subset(results_subset, results_subset$TDG09_sites == 0)

length(with_sites$goTerm_frac)  # 921
length(without_sites$goTerm_frac)  # 2856

par(mfrow=c(1,2))
hist(with_sites$goTerm_frac)
hist(without_sites$goTerm_frac)

ks.test(with_sites$goTerm_frac,without_sites$goTerm_frac)  # p-value of 0.02518

hist(with_sites$goTerm_frac[with_sites$goTerm_frac != 0], ylim=c(0,3.0), main=">= 1 sites found by TDG", xlab="GO-term fraction", freq=F)  # 125
hist(without_sites$goTerm_frac[without_sites$goTerm_frac != 0], main="0 sites found by TDG", xlab="GO-term fraction",freq=F)  # 449
