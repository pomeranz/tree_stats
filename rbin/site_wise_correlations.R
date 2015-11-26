# paml, slr and fubar postmean and variance extraction from data. 

# Packages
library(yaml)
library(iterators)
library(stringi)
library(stringr)
library(base)
library(xtable)
library(ape)

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
  
  file <- readLines(infile)
  
  # Error 3 - unfinished output --> job was killed half way through
  if (!("(BEB)" %in% unlist(str_split(file, pattern= " ")))) {
    run = FALSE
    return(list("Error: unfinished job", run))
  }
  
  
  # locating correct part of the file
  beb <- grep("(BEB)", file)
  end_of_file <- length(file)
  sub_file <- file[beb:end_of_file]
  site_wise <- 4  # checked manually after how many lines the site-wise values appear
  end_of_subfile <- length(sub_file)
  positive <- grep("Positively", sub_file)+4
  end_site_wise <- positive-6  # checked manually where the end of the values are 
                                # from the the start of the positively selected sites
  site_wise_values <- sub_file[site_wise:end_site_wise]
  
  # setting up loop
  index = 1
  df <- data.frame(matrix(ncol=3))
  site_wise_df <- data.frame(matrix(ncol=3))
  
  # not working yet. 
  for (i in site_wise_values) {
    # fix the few ones whose vector length is one less
    temp = unlist(str_split(i, pattern=" ")) 
    if ("(11)" %in% temp) {
      weird = which(temp == "(11)")
      temp = c(temp[1:weird], "filler", temp[(weird+1):length(temp)])
      i = paste(temp, collapse=" ")
    }
  
    if (index %in% seq(1,9)) {
      site <- unlist(str_split(i, pattern=" "))[4]
      value <- unlist(str_split(i, pattern=" "))[22]
      variance <- unlist(str_split(i, pattern=" "))[25]
    } else if (index %in% seq(10,99)) {
      site <- unlist(str_split(i, pattern=" "))[3]
      value <- unlist(str_split(i, pattern=" "))[21]
      variance <- unlist(str_split(i, pattern=" "))[24]
    } else if (index %in% seq(100,999)) {
      site <- unlist(str_split(i, pattern=" "))[2]
      value <- unlist(str_split(i, pattern=" "))[20]
      variance <- unlist(str_split(i, pattern=" "))[23]
    } else {
      site <- unlist(str_split(i, pattern=" "))[1]
      value <- unlist(str_split(i, pattern=" "))[19]
      variance <- unlist(str_split(i, pattern=" "))[22]
    }
    
    site_wise_df[index,] <- c(site,value,variance)
    index = index + 1
  }
  
  index = 1
  
  for (i in sub_file[positive:end_of_subfile]) {  # this is for the positiveley selected sites
  
    if (length(unlist(str_split(i, pattern= " "))) > 30) {
      break
    } else {
      site <- as.numeric(substr(i, 1,7))
      prob <- as.numeric(substr(i, 14,19))
      omega <- substr(i,29,43)
      # maybe make the probability check here already
      if (prob >= 0.95) {
        df[index,] <- c(site,prob,omega)
      }
      index = index + 1
    }
  }
  
  
  
  colnames(df) <- c("site", "Prob(w>1)", "mean w")
  colnames(site_wise_df) <- c("site", "value", "variance")
  
  run = TRUE
  
  # get the likelihoods to calculate wether it is even worth looking at the alignment
  
  lnLs <- grep("lnL", file)
  model_eight <- lnLs[2]
  model_seven <- lnLs[1]
  deltalnL <- as.numeric(substr(file[model_eight], 6, 20)) - as.numeric(substr(file[model_seven], 6, 20))
  p_val <- pchisq(2*deltalnL, df=2, lower.tail=FALSE)
  lnL_info <- c(infile, as.numeric(substr(file[model_seven], 6, 20)), as.numeric(substr(file[model_eight], 6, 20)), deltalnL, p_val)
  
  
  return(list(df, run, nrow(df), lnL_info, site_wise_df))
  # function end  
}

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

########################################## FUBAR results extraction function ##################################################

fubar_results <- function(infile) {
  
  if (file.exists(infile) == FALSE) {
    run = FALSE
    return(list("No data available", run, 0))
  }
  
  fubar_result <- read.csv(infile)
  
  run = TRUE
  
  # extract the significant hits some other time.
  
  # calculate dN/dS (beta/alpha)
  
  omega = apply(fubar_result, 1, FUN=function(x) {x[3]/x[2]})
  
  fubar_result$omega <- omega
  
  
  return(list(fubar_result,run))
# end of function  
}

########################################## Plot saving fucntion ##################################################

plot_create <- function(name, subdir, basename, points1, points2) {
  
  if (name == "PAML_SLR") {
    setwd(paste("/home/gideon/Documents/mphil_internship/results_out/correlation_analysis/PAML_SLR/", subdir, "/",sep = ""))
    pdf(basename)
    plot(points1,points2, xlab="PAML_omega", ylab="SLR_omega", main = paste("Correlation of omega values PAML_SLR of file", basename, sep=" "))
    abline(a=0, b=1, col="red")
    if(any(points1 >= 1)) {
      rect(1,0,max(points1, na.rm=T), max(points2, na.rm=T), col=rgb(1,0,0,alpha=0.1), border=NA)
    }
    if (any(points2 >= 1)) {
      rect(0,1,max(points1, na.rm=T), max(points2, na.rm=T), col=rgb(1,0,0,alpha=0.1), border=NA)
    }
    dev.off()
  } else if (name == "PAML_FUBAR") {
    setwd(paste("/home/gideon/Documents/mphil_internship/results_out/correlation_analysis/PAML_FUBAR/", subdir, sep = ""))
    pdf(basename)
    plot(points1,points2, xlab="PAML_omega", ylab="FUBAR_omega", main = paste("Correlation of omega values PAML_FUBAR of file", basename, sep=" "))
    abline(a=0, b=1, col="red")
    if (any(points1 >= 1)) {
      rect(1,0,max(points1, na.rm=T), max(points2, na.rm=T), col=rgb(1,0,0,alpha=0.1), border=NA)
    }
    if (any(points2 >= 1)) {
      rect(0,1,max(points1, na.rm=T), max(points2, na.rm=T), col=rgb(1,0,0,alpha=0.1), border=NA)
    }
    dev.off()
  } else {
    setwd(paste("/home/gideon/Documents/mphil_internship/results_out/correlation_analysis/SLR_FUBAR/", subdir, sep = ""))
    pdf(basename)
    plot(points1,points2, xlab="SLR_omega", ylab="FUBAR_omega", main = paste("Correlation of omega values SLR_FUBAR of file", basename, sep=" "))
    abline(a=0, b=1, col="red")
    if (any(points1 >= 1)) {
      rect(1,0,max(points1, na.rm=T), max(points2, na.rm=T), col=rgb(1,0,0,alpha=0.1), border=NA)
    }
    if (any(points2 >= 1)) {
      rect(0,1,max(points1, na.rm=T), max(points2, na.rm=T), col=rgb(1,0,0,alpha=0.1), border=NA)
    }
    dev.off()
  }
} 


# make it split into 4 parts above 0 and then only one category for lower than 0 
plot_save <- function(name,corr,basename, points1, points2) {
  est <- corr$estimate
  if (name == "PAML_SLR") {
    if (est >= 0.75) {  # positive
      plot_create(name="PAML_SLR",subdir="positive", basename, points1, points2)
    } else if (est >= 0.50) {  # good 
      plot_create(name="PAML_SLR",subdir="good", basename, points1, points2)
    } else if (est >= 0.25) {  # ok
      plot_create(name="PAML_SLR",subdir="ok", basename, points1, points2)
    } else if (est >= 0.0) {  # neutral
      plot_create(name="PAML_SLR",subdir="neutral", basename, points1, points2)
    } else {  # bad
      plot_create(name="PAML_SLR",subdir="negative", basename, points1, points2)
    }
    # end of PAML_SLR  
  } else if (name == "PAML_FUBAR") {
    if (est >= 0.75) {  # positive
      plot_create(name="PAML_FUBAR",subdir="positive", basename, points1, points2)
    } else if (est >= 0.50) {  # good 
      plot_create(name="PAML_FUBAR",subdir="good", basename, points1, points2)
    } else if (est >= 0.25) {  # ok
      plot_create(name="PAML_FUBAR",subdir="ok", basename, points1, points2)
    } else if (est >= 0.0) {  # neutral
      plot_create(name="PAML_FUBAR",subdir="neutral", basename, points1, points2)
    } else {  # bad
      plot_create(name="PAML_FUBAR",subdir="negative", basename, points1, points2)
    }
    # end of PAML_FUBAR
  } else {
    if (est >= 0.75) {  # positive
      plot_create(name="SLR_FUBAR",subdir="positive", basename, points1, points2)
    } else if (est >= 0.50) {  # good 
      plot_create(name="SLR_FUBAR",subdir="good", basename, points1, points2)
    } else if (est >= 0.25) {  # ok
      plot_create(name="SLR_FUBAR",subdir="ok", basename, points1, points2)
    } else if (est >= 0.0) {  # neutral
      plot_create(name="SLR_FUBAR",subdir="neutral", basename, points1, points2)
    } else {  # bad
      plot_create(name="SLR_FUBAR",subdir="negative", basename, points1, points2)
    }
    # end of SLR_FUBAR
  }
# end of function
}


########################################## LOOP ##################################################

# The loop extractes the site-wise values for all 3 methods and produces 3 plots at each iteration. PAML-SLR, PAML-FUBAR, SLR-FUBAR
# based on the correlation coefficient I put them into a few categories: 
# positive(0.75 - 1), 
# good(0.50-0.74), 
# ok(0.25-0.49), 
# neutral(0.0-0.24), 
# negative(-0.1- -1)


# Setup
positive <- c()
ok <- c()
neutral <- c()
bad <- c()
negative <- c()
# table to hold all the info as well
corr_table <- data.frame(matrix(ncol=7))
colnames(corr_table) <- c("file", "PAML_SLR_corr", "PAML_SLR_p_val", "PAML_FUBAR_corr", "PAML_FUBAR_p_val", "SLR_FUBAR_corr", "SLR_FUBAR_p_val")

# create folders 3 x 5
setwd("/home/gideon/Documents/mphil_internship/results_out")
dir.create(path="/home/gideon/Documents/mphil_internship/results_out/correlation_analysis/")
savedir = "/home/gideon/Documents/mphil_internship/results_out/correlation_analysis/"

# categories
categ <- c("PAML_SLR", "PAML_FUBAR", "SLR_FUBAR")
# subcategories
subcateg <- c("positive", "good", "ok", "neutral", "negative")

# create the directories
for (i in categ) {
  dir.create(path=paste(savedir,i, sep=""))
  for (j in subcateg) {
    dir.create(path=paste(savedir,i, "/",j, sep=""))
  }
}


# index
index = 1

setwd("/home/gideon/Documents/mphil_internship/")
inroot = getwd()

for (infile in Sys.glob(file.path(inroot, "prank_out", "*/*_prank.best.fas"))) {
  
  print(infile)
  
  # create basenames to acces the other files.
  basename = head(unlist(strsplit(tail(unlist(strsplit(infile, split="/")), n=1), split=".", fixed=TRUE)), n=1)
  basename = paste(unlist(strsplit(basename, split="_"))[1], unlist(strsplit(basename, split="_"))[2], sep="_")
  prefix = substr(unlist(strsplit(basename, split="_", fixed=TRUE))[1], 1, 2)
  
  aln_name <- paste(basename, ".fa", sep="")
  
  corr_table[index,1] <- aln_name
  
  paml_file <- paste("/home/gideon/Documents/mphil_internship/paml_out/out/",prefix, "/", basename, "/", "rst", sep="")
  paml_res_plus <- paml_results(paml_file)
  paml_res <- paml_res_plus
  
  slr_file <- paste("/home/gideon/Documents/mphil_internship/slr_out/Eutheria/",prefix, "/", basename, "_matched.res", sep="")
  slr_res_plus <- slr_results(slr_file)
  
  fubar_file <- paste("/home/gideon/Documents/mphil_internship/fubar_prep_out/Eutheria/",prefix, "/", basename, ".nwk.fubar.csv", sep="")
  fubar_res_plus <- fubar_results(fubar_file)
  
  
  if (paml_res_plus[[2]] == TRUE & slr_res_plus[[2]] == TRUE & fubar_res_plus[[2]] == TRUE) {
    paml_points <- as.numeric(paml_res_plus[[5]]$value)
    slr_points <- as.numeric(slr_res_plus[[4]]$omega)
    fubar_points <- as.numeric(fubar_res_plus[[1]]$omega)
    
    
    # remove the points that SLR marked as single char
    to_remove <- slr_res_plus[[5]]
    if (length(to_remove) != 0) {
      slr_points <- slr_points[-to_remove]
      paml_points <- paml_points[-to_remove]
      fubar_points <- fubar_points[-to_remove]
    }
    
    # do correlations and plots and put into correct file path
    # PAML_SLR
    
    p_s_corr <- cor.test(paml_points,slr_points)
    corr_table[index,2] <- p_s_corr$estimate
    corr_table[index,3] <- p_s_corr$p.value
    
    plot_save(name="PAML_SLR", p_s_corr, basename, paml_points, slr_points)
    
    # PAML_FUBAR
    p_s_corr <- cor.test(paml_points,fubar_points)
    corr_table[index,4] <- p_s_corr$estimate
    corr_table[index,5] <- p_s_corr$p.value
    
    plot_save(name="PAML_FUBAR", p_s_corr, basename, paml_points, fubar_points)
    
    # SLR_FUBAR
    p_s_corr <- cor.test(slr_points,fubar_points)
    corr_table[index,6] <- p_s_corr$estimate
    corr_table[index,7] <- p_s_corr$p.value
    
    plot_save(name="SLR_FUBAR", p_s_corr, basename, slr_points, fubar_points)
       
  } else {
    corr_table[index,2:7] <- NA
  }
  
  index = index + 1
# end of loop  
}

################################################# Analysis ##############################################################

setwd("/home/gideon/Documents/mphil_internship/results_out/correlation_analysis/")
write.csv(corr_table, file="corr_table.csv")
