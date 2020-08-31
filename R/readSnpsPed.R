#Conversion to single column per locus from plink file via LEA function
# readSnpsPed <- function (infile, inames, lnames){
#   if (missing(inames))
#     stop("missing 'inames' file")
#   if (missing(lnames))
#     stop("missing 'lnames' file")
#
#   LEA::ped2lfmm(infile, output.file=tempfile, force=TRUE)
#   snpobj <- read.table(tempfile)
# 
#   
#   rownames(snpobj) <- inames
#   colnames(snpobj) <- lnames
#   
#   return(snpobj)
# }

#Two columns per locus, direct import from plink files
readSnpsPed <- function (pedfile, mapfile){
  if (missing(pedfile))
    stop("missing plink '.ped' file")
  if (missing(mapfile))
    stop("missing plink '.map' file")
  
  #Isolate SNP matrix from .ped file and and append sample IDs as row names
  snpobj <- read.table(pedfile, row.names=2, na.strings = c("0", "-9"), stringsAsFactors = FALSE)
  snpobj <- snpobj[-(1:5)]
  
  #Extract locus IDs from .map file, create column names and append to SNP matrix
  lnames <- read.table(mapfile)[,2]
  locnames1 <- paste0(lnames, ".1")
  locnames2 <- paste0(lnames, ".2")
  locnames3 <- c(rbind(locnames1,locnames2))
  colnames(snpobj) <- locnames3
  
  if (any(is.na(snpobj))){
    hist(colSums(is.na(as.matrix(snpobj))), main="Missing genotypes per SNP", xlab="No. missing genotpyes")
    Sys.sleep(1)
    impute <- readline("This dataset contains NAs. Please select an imputation method (type a number): 0. none, 1. mode ")
    
    if (impute == "1"){
      
      #Function for finding the major allele for each column
      Mode <- function(x, na.rm = FALSE) {
        if(na.rm){
          x = x[!is.na(x)]
        }
    
        ux <- unique(x)
        return(ux[which.max(tabulate(match(x, ux)))])
      }
  
      for(i in 1:ncol(snpobj)){
        snpobj[i] <- na.replace(snpobj[i], Mode(snpobj[i], na.rm=TRUE)) #replace NAs with the major allele at a given locus
        snpobj[i] <- replace(snpobj[i], snpobj[i] == Mode(snpobj[i], na.rm=TRUE), "0") #recode the major allele as "0"
        snpobj[i] <- replace(snpobj[i], (snpobj[i] == "A" | snpobj[i] == "G" | snpobj[i] == "T" | snpobj[i] == "C"), "1") #recode the minor allele as "1"
      }
    }
  }  else{
    
    for(i in 1:ncol(snpobj)){
    snpobj[i] <- replace(snpobj[i], snpobj[i] == Mode(snpobj[i], na.rm=TRUE), "0") #recode the major allele as "0"
    snpobj[i] <- replace(snpobj[i], (snpobj[i] == "A" | snpobj[i] == "G" | snpobj[i] == "T" | snpobj[i] == "C"), "1") #recode the minor allele as "1"
    }
  }
  
  return(data.matrix(snpobj, rownames.force = TRUE))

}