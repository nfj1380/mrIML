#'Calculate binary outcome deviance residuals using probability predictions
#'
#'
#'This function uses the predictions from a model applied to binary outcomes to
#'calculate deviance residuals. Predictions must be on the probability scale.
#'
#'@param yhat A \code{vector} of predictions on the probability scale.
#'@param truth A \code{vector} of true binary observations in the same order as \code{yhat}
#'@export
#'
devianceResids = function(yhat, truth){
  # Check for appropriate data structures
  if(any(class(truth) %in% 'tbl')){
    truth <- dplyr::pull(truth, 1)
  }
  
  if(any(class(yhat) %in% 'tbl')){
    yhat <- dplyr::pull(yhat, 1)
  }
  
  if(any(is.na(yhat))){
    stop('NAs detected in yhat')
  }
  
  if(any(is.na(truth))){
    stop('NAs detected in truth')
  }
  
  is.binary <- function(v) {
    x <- unique(v)
    length(x) - sum(is.na(x)) == 2L && all(sort(x[1:2]) == 0:1)
  }
  
  if(!is.binary(truth)){
    stop('truth must be a binary vector (1s and 0s only)')
  }
  
  # To calculate deviance residuals consistently for binomial models, we need a simple formula
  # truth equals binary observations; yhat equals model probability predictions
  resid <- ifelse(truth == 1, 
                  sqrt((-2 * log(yhat))), 
                  -1 * (sqrt((-2*log(1 - yhat)))))
  resid <- as.vector(resid)
  return(resid)
}
