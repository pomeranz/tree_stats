# fubar results extraction function
# extra function to find out extremes of alpha. 

########################################## FUBAR results extraction function ##################################################

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


# testing
setwd("/home/gideon/Documents/mphil_internship/")
inroot = getwd()

# setup
fubar_temp <- data.frame(matrix(ncol=5))
colnames(fubar_temp) <- c("name", "sites", "sd(alpha)", "sd(beta)", "sd(omega)")

index = 1

for (infile in Sys.glob(file.path(inroot, "prank_out", "*/*_prank.best.fas"))) {
  
  print(infile)
  
  basename = head(unlist(strsplit(tail(unlist(strsplit(infile, split="/")), n=1), split=".", fixed=TRUE)), n=1)
  basename = paste(unlist(strsplit(basename, split="_"))[1], unlist(strsplit(basename, split="_"))[2], sep="_")
  prefix = substr(unlist(strsplit(basename, split="_", fixed=TRUE))[1], 1, 2)
  
  fubar_temp[index,1] <- paste(basename, sep="")
  
  fubar_file <- paste("/home/gideon/Documents/mphil_internship/fubar_prep_out/Eutheria/",prefix, "/", basename, ".nwk.fubar.csv", sep="")
  
  fubar_res_plus <- fubar_results(fubar_file)
  fubar_sites <- fubar_res_plus[[3]]
  
  if (fubar_res_plus[[2]] == FALSE) {
    fubar_temp[index,2:5] <- NA
    
  } else {
    fubar_res <- fubar_res_plus[[4]]
    fubar_temp[index,2] <- fubar_sites
    fubar_temp[index,3] <- sd(fubar_res$alpha)
    fubar_temp[index,4] <- sd(fubar_res$beta)
    fubar_temp[index,5] <- sd(fubar_res$omega)
  }
  
  
  index = index + 1
}

# Analysis

mean(fubar_temp$X2, na.rm=T)
max(fubar_temp$X2, na.rm=T)

fubar_temp[which(fubar_temp$X2 != 0),]  # 2608

which(is.na(fubar_temp[,2]))  # 1, wow just one
