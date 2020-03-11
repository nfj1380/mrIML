#' Utility function to create multi-fold cross-validation with each fold balanced 
#'
#' This function allows your dataframe to be be partitioned for multi-fold cross-validation with a balanced number of classes (e.g.same number of presences absences). Called as part of the mcaret function.
#' @param y is a data.frame 
#' @param k is the number of partiions (default is 10) 
#' @param times is the number of times the data is partioned (default is 5) 
#' @keywords cross-validation
#' @export

CreateBalancedMultiFolds <- function (y, k = 10, times = 5) 
{
  prettyNums <- paste("Rep", gsub(" ", "0", format(1:times)), 
                      sep = "")
  for (i in 1:times) {
    tmp <- CreateBalancedFolds(y, k = k, list = TRUE, returnTrain = TRUE)
    names(tmp) <- paste("Fold", gsub(" ", "0", format(seq(along = tmp))), 
                        ".", prettyNums[i], sep = "")
    out <- if (i == 1) 
      tmp
    else c(out, tmp)
  }
  out
}
