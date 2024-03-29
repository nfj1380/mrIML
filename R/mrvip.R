#' Wrapper to estimate model-agnostic variable importance for multi-response models.
#'
#' @param yhats A \code{list} is the list generated by mrIMLpredicts
#' @param X  A \code{dataframe} is the feature data 
#' @param Y A \code{dataframe} is the response data
#' @return A dataframe containing variable importance scores for each response variable
#' 
#' @details Calculates variable importance using based on partial dependencies but could do permutations as well (more memory intensive).
#' Key input is the object created by MrIML. Can be plotted with the plot_vi function.
#' @examples 
#' VI <- mrVip(yhats, X=X, Y=Y)
#' groupCov <- c(rep ("Host_characteristics", 1),rep("Urbanisation", 3), rep("Vegetation", 2), rep("Urbanisation",1), rep("Spatial", 2), 
#' rep('Host_relatedness', 6),rep ("Host_characteristics", 1),rep("Vegetation", 2), rep("Urbanisation",1))  
#' plot_vi(VI=VI,  X=X,Y=Y, modelPerf=ModelPerf, groupCov, cutoff= 0.5)
#' 
#' @export 
#' 
mrVip <- function (yhats, X, Y) { 
  
  n_response <- length(yhats)
  n_features <- sort(names(X))
  
  modList <- yhats %>% purrr::map(pluck('mod1_k')) # extracts model
  
  trainList <- yhats %>% purrr::map(pluck('data_train'))
  
  imp_list <- lapply(seq_along(Y), function(i) {
    
    resp <- names(Y)[i]
    
    imp <- modList[[i]] %>%  
      extract_fit_parsnip()
    
    imp$fit$coefficients[is.na(imp$fit$coefficients)] <- 0 # some algorithms will bring back NA coefficients
    
    #hstats way
    # impVIcalc <- hstats(imp, X = X)
    # impVI <- pd_importance(impVIcalc)
   
    impVI <- vi(imp, method='firm',feature_names= names(X), train=trainList[[i]][-1] )
    #impVI <- vi(imp, num_features = length(X)) 
  
    impD <- impVI %>% 
      arrange(Importance) %>% #changed from Variable
      pluck("Importance")
    
    missing <- imp$fit$coefficients[is.na(imp$fit$coefficients)] 
    missing[is.na(missing)] <- 0
    impDcombined <- c(impD, missing)
    names(impDcombined) <- n_features
    
    imp_df <- as.data.frame(impDcombined)
    names(imp_df) <- resp
    
    return(imp_df)
  })
  
  ImpGlobal <- bind_cols(!!!imp_list)
  
  return(ImpGlobal)
}
