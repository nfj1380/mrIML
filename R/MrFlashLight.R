mrFlashlight <- function(yhats, X, Y, response = "multi", index = 1, model = "regression") {
  
  if (model == "classification") {
    
    
    
    # Prediction function
    
    pred_fun <- function(m, dat) {
      
      predict(
        
        m,
        
        data.matrix(dat[, colnames(Y), drop = FALSE]),
        
        type = "prob"
        
      )[[".pred_1"]]
      
    }
    
    
    
    # List of metrics
    
    metrics = list(
      
      logloss = MetricsWeighted::logLoss,
      
      `ROC AUC` = MetricsWeighted::AUC,
      
      `% Dev Red` = MetricsWeighted::r_squared_bernoulli
      
    )
    
    
    
    # Create explainer
    
    if (response == "single") {
      
      mfl <- flashlight(
        
        model = yhats[[index]]$mod1_k, #change the index to focus on other SNPs
        
        label = colnames(X)[index],
        
        data = cbind(Y, X),
        
        y = colnames(X)[index],
        
        predict_function = pred_fun,
        
        metrics = metrics
        
      )
      
    } else if(response == "multi") {
      
      models <- lapply(yhats, `[[`, "mod1_k")
      
      fl_list <- vector("list", length(yhats))
      
      for (i in seq_along(fl_list)) {
        
        fl_list[[i]] <- flashlight(
          
          model = yhats[[i]]$mod1_k,
          
          label = colnames(X)[i],
          
          y = colnames(X)[i]
          
        )
        
      }
      
      mfl <- multiflashlight(
        
        fl_list,
        
        data = cbind(Y, X),
        
        predict_function = pred_fun,
        
        metrics = metrics
        
      )
      
    }
    
  }
  
  
  
  if (model == "regression") { 
    
    
    
    # Prediction function
    
    pred_fun <- function(m, dat) {
      
      predict(m, data.matrix(dat[, colnames(Y), drop = FALSE]))[[".pred"]]
      
    }
    
    
    
    # List of metrics
    
    metrics = list( 
      
      rmse = MetricsWeighted::rmse,
      
      `R-squared` = MetricsWeighted::r_squared
      
    )
    
    
    
    # Create explainer
    
    if (response == "single") {
      
      mfl <- flashlight(
        
        model = yhats[[index]]$mod1_k, #change the index to focus on other SNPs
        
        label = colnames(X)[index],
        
        data = cbind(Y, X),
        
        y = colnames(X)[index],
        
        predict_function = pred_fun,    #predict_function = pred_fun,
        
        metrics = metrics
        
      )
      
    } else if (response == "multi") {
      
      models <- lapply(yhats, `[[`, "mod1_k")
      
      fl_list <- vector("list", length(yhats))
      
      for (i in seq_along(fl_list)) {
        
        fl_list[[i]] <- flashlight(
          
          model = yhats[[i]]$mod1_k,
          
          label = colnames(X)[i],
          
          y = colnames(X)[i]
          
        )
        
      }
      
      mfl <- multiflashlight(
        
        fl_list,
        
        data = cbind(Y, X),
        
        predict_function = pred_fun,
        
        metrics = metrics
        
      )
      
    }
    
  }
  
  return(mfl)
}