readSnpsPed <- function (infile, inames, lnames){
  if (missing(inames))
    stop("missing 'inames' file")
  if (missing(lnames))
    stop("missing 'lnames' file")
  tempfile <- "temp.lfmm"
  LEA::ped2lfmm(infile, output.file=tempfile, force=TRUE)
  snpobj <- read.table(tempfile)
  file.remove(tempfile)
  
  rownames(snpobj) <- inames
  colnames(snpobj) <- lnames
  
  return(snpobj)
}