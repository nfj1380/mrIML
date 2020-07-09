#' Utility function to create cross-validation with each fold balanced 
#'
#' This function allows your dataframe to be be partitioned for multi-fold cross-validation with a balanced number of classes (e.g.same number of presences absences). Called as part of the mcaret function.
#' @param y is a data.frame #add more Gustavo?
#' @param k is the number of partiions (default is 10) 
#' @keywords cross-validation
#' @export

CreateBalancedFolds <- function (y, k = 10, list = TRUE, returnTrain = FALSE) 
{
  if (is.numeric(y)) {
    cuts <- floor(length(y)/k)
    if (cuts < 2) 
      cuts <- 2
    if (cuts > 5) 
      cuts <- 5
    y <- cut(y, unique(quantile(y, probs = seq(0, 1, length = cuts))), 
             include.lowest = TRUE)
  }
  
  if (k < length(y)) {
    y <- factor(as.character(y))
    numInClass <- table(y)
    foldVector <- vector(mode = "integer", length(y))
    
    # balance classes to sample examples in CV
    numInClass[1:length(numInClass)]<-min(numInClass)
    
    for (i in 1:length(numInClass)) {
      seqVector <- rep(1:k, numInClass[i]%/%k)
      if (numInClass[i]%%k > 0) 
        seqVector <- c(seqVector, sample(1:k, numInClass[i]%%k))
      foldVector[sample(which(y == dimnames(numInClass)$y[i]),numInClass[1])] <- sample(seqVector)
    }
  }
  else foldVector <- seq(along = y)
  
  if (list) {
    out <- split(seq(along = y), foldVector)
    # remove `0` element: examples not in balanced data set
    out <- out[-1]
    # save ids for examples of the balanced dataset
    idxBalanced <- unique(unlist(out))
    names(out) <- paste("Fold", gsub(" ", "0", format(seq(along = out))), sep = "")
    # return indexes for training data, removing the testing examples (out$) from the balanced 
    # data examples (idxBalance), for each fold
    if (returnTrain) 
      out <- lapply(out, function(data, y) y[setdiff(idxBalanced,data)], y = seq(along = y))
  }
  else out <- foldVector
  out
}


