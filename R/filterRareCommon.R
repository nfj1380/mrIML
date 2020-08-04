#' Filter rare or common OTUs or SNPs or species from your response data
#'
#' This function allows you to remove OTUs or SNPs or species from your response data as a preprocessing step. Suitable when the response is a binary outcome.
#' @param X is a data.frame with rows as sites or individuals or populations and columns as loci or species OTUs.
#' @param lower is the lower threshold value  in which response varialkes are removed from the data.frame.
#' @param higher is the upper threshold value  in which response varialkes are removed from the data.frame.
#' @keywords filter
#' @export

filterRareCommon <- function(X, lower=lower, higher=higher){
  n = ncol(X)
  r = nrow(X)
  Xt <- as.data.frame(t(X))
  Xt$Xsum <- rowSums(Xt[1:n,] )
  FilterNoCommon <-subset(Xt, Xsum < (r*0.6)) 
  FilterNoRare <-subset(FilterNoCommon , Xsum > (r*0.4))
  #remove 'new' sort column
  FilterNoRare$Xsum <- NULL
  X <- as.data.frame(t(FilterNoRare))
  X
}