# simulation results

# Packages 
library(yaml)
library(iterators)
library(stringi)
library(stringr)
library(base)
library(xtable)
library(ape)


########################################## SLR results extraction function ##################################################

slr_results <- function(infile) {
  
  if (file.exists(infile) == FALSE) {
    run = FALSE
    return(list("No data available", run, 0))
  }
  
  
  slr_result <- readLines(infile)
  slr_result <- slr_result[-1]
  
  single_char <- grep("Single char", slr_result)
  all_gaps <- grep("All gaps", slr_result)
  
  
  cn <- c("Site", "Neutral", "Optimal", "omega", "lower", "upper", "LRT_Stat", "Pval", "Adj.Pval", "Q-value")
  
  slr_result <- read.table(infile, fill = TRUE, row.names=NULL, header = FALSE)
  
  if (any(slr_result[,8] != "NA")) {
    slr_result <- slr_result[which(slr_result[,8] != "NA"),]
  }
  
  
  if (any(slr_result[,1] == "char")) {
    slr_result <- slr_result[-which(slr_result[,1] == "char"),]
  }
  
  if (any(slr_result[,1] == "gaps")) {
    slr_result <- slr_result[-which(slr_result[,1] == "gaps"),]
  }
  
  if (any(slr_result[,1] == "Constant")) {
    slr_result <- slr_result[-which(slr_result[,1] == "Constant"),]
  }
  
  if (any(slr_result[,1] == "Synonymous")) {
    slr_result <- slr_result[-which(slr_result[,1] == "Synonymous"),]
  }
  
  if (any(slr_result[,1] == "Single")) {
    slr_result <- slr_result[-which(slr_result[,1] == "Single"),]
  }
  
  if (any(slr_result[,1] == "!")) {
    slr_result <- slr_result[-which(slr_result[,1] == "!"),]
  }
  
  if (any(slr_result[,1] == "All")) {
    slr_result <- slr_result[-which(slr_result[,1] == "All"),]
  }
  
  rownames(slr_result) <- NULL
  slr_result <- slr_result[,1:10]
  
  colnames(slr_result) <- cn
  
  #slr_result$omega <- gsub(" ", "", substr(slr_results$omega, 1, nchar(slr_results$omega)-1), fixed=TRUE)
  #slr_result$omega <- as.numeric(slr_results$omega)
  #slr_result$Adj.Pval <- gsub(" ", "", substr(slr_results$Adj.Pval, 1, nchar(slr_results$Adj.Pval)-1), fixed=TRUE)
  #slr_result$Adj.Pval <- as.numeric(slr_results$Adj.Pval)
  
  #colnames(slr_result) <- cn
  
  
  signif_vals <- which(slr_result$Adj.Pval <= 0.05 & slr_result$omega > 1.0)  # Which hits
  no_signifs <- length(signif_vals)
  
  output <- slr_result[signif_vals,]
  
  final_out <- output[,-c(1,2,4,5)]
  #rownames(final_out) <- gsub(" ", "", substr(rownames(final_out), 1, nchar(rownames(final_out))-1), fixed=TRUE)
  
  run = TRUE
  
  return(list(final_out,run,length(rownames(slr_result)), slr_result, c(single_char,all_gaps)))
  # end of function
}
################################### FUBAR EXTRACTION FUNCTION ###############################################

fubar_results <- function(infile) {
  
  if (file.exists(infile) == FALSE) {
    run = FALSE
    return(list("No data available", run, 0))
  }
  
  fubar_result <- read.csv(infile)
  
  run = TRUE
  
  # calculate dN/dS (beta/alpha)
  
  omega = apply(fubar_result, 1, FUN=function(x) {x[3]/x[2]})
  
  fubar_result$omega <- omega
  
  # extract the significant hits some other time.
  signif <- fubar_result[which(fubar_result[,5] >= 0.95 & fubar_result$omega >= 1),]
  no_signifs <- nrow(signif)
  
  
  return(list(signif,run, no_signifs, fubar_result))
  # end of function  
}


############################## LOOP ##################################

# for fubar plot the dS values. 
# run1 had dS values of 0.5 (10% of time) and 1.0 (90% of time)
# plot dS for each alignment and create a histogram. for each tree seperatly

# for SLR plot omega values and note in legend how many sites under positive selection. 

species_numbers <- c("6species", "12species", "17species", "44species")
sizes <- c("Big", "Medium", "Small")

### MAKE SURE TO CHECK THAT ALL DIRECTORIES ARE IN PLACE. MIGHT NEED MANUAL REFINEMENT

# create the results and run directories yourself
savedir_fubar <- "/home/gideon/Documents/mphil_internship/simulations/results/run1/fubar"
dir.create(savedir_fubar)
savedir_slr <- "/home/gideon/Documents/mphil_internship/simulations/results/run1/slr"
dir.create(savedir_slr)

savedirs <- c(savedir_fubar, savedir_slr)

# create the directories
for (savedir in savedirs) {
  for (i in species_numbers) {
    dir.create(path=file.path(savedir,i))
    for (j in sizes) {
      dir.create(path=file.path(savedir,i,j))
    }
  }
}


for (species in species_numbers) {
  
  print(species)
  
  for (size in sizes) {
    
    print(size)
    
    index = 1
    
    
    for (infile in Sys.glob(file.path("/home/gideon/Documents/mphil_internship/simulations/fubar_prep_out1", species, size, "*", "*.csv"))) {
      
      fubar_res_plus <- fubar_results(infile)
      
      alpha <- fubar_res_plus[[4]]$alpha
      
      setwd(file.path(savedir_fubar, species, size))
      png(index)
      hist(alpha,breaks=100, xlab = "alpha")
      dev.off()
      
      # create basenames to acces the other files.
      basename = head(unlist(strsplit(tail(unlist(strsplit(infile, split="/")), n=1), split=".", fixed=TRUE)), n=1)
      prefix = basename
      
      
      slr_file <- file.path("/home/gideon/Documents/mphil_internship/simulations/slr_prep_out1",species, size, prefix, paste(basename, "_matched.res", sep=""))
      
      slr_res_plus <- slr_results(slr_file)
      
      omega <- slr_res_plus[[4]]$omega
      no_sites <- nrow(slr_res_plus[[1]])
      
      setwd(file.path(savedir_slr, species, size))
      png(index)
      hist(omega,breaks=100, xlab = "omega")
      legend("topright", legend=paste("pos_sel sites:", no_sites))
      dev.off()
      
      
      index = index + 1
      
    }
    
  # end of size loop  
  }
  
# end of species loop 
}
