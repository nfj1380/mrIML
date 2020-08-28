#Conversion to single column per locus from plink file via LEA function
# readSnpsPed <- function (infile, inames, lnames){
#   if (missing(inames))
#     stop("missing 'inames' file")
#   if (missing(lnames))
#     stop("missing 'lnames' file")
#   tempfile <- "temp.lfmm"
#   LEA::ped2lfmm(infile, output.file=tempfile, force=TRUE)
#   snpobj <- read.table(tempfile)
#   file.remove(tempfile)
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
  snpobj <- read.table(pedfile, row.names=2)
  snpobj <- snpobj[-(1:5)]
  
  #Extract locus IDs from .map file, create column names and append to SNP matrix
  lnames <- read.table(mapfile)[,2]
  locnames1 <- paste0(lnames, ".1")
  locnames2 <- paste0(lnames, ".2")
  locnames3 <- c(rbind(locnames1,locnames2))
  colnames(snpobj) <- locnames3
  
  return(snpobj)
}