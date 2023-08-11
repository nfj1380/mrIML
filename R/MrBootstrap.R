

mrBootstrap <- function(yhats, num_bootstrap = 10, Y=Y,  alpha = 0.05, ice = FALSE) {
  
  n_response <- length(yhats)
  
  pb <- txtProgressBar(min = 0, max = n_response, style = 3) 
  
  internal_fit_function <- function(k) {
    
    setTxtProgressBar(pb, k) #progressbar marker
    
    features <- colnames(yhats[[k]]$data)[-1]
    
    n <- nrow(yhats[[k]]$data)
    
    pd_raw <- vector("list", num_bootstrap)  # Initialize pd_raw as a list
    
    for (i in 1:num_bootstrap) {
      # Generate bootstrap sample
      bootstrap_sample <- yhats[[k]]$data[sample(1:n, replace = TRUE), ] ###
      
      # Extract the workflow from the best fit
      wflow <- yhats[[k]]$last_mod_fit %>% extract_workflow()
      
      # Add the bootstrap data to the workflow
      wflow$data <- bootstrap_sample
      
      # Fit the model using the bootstrap sample
      model_fit <- fit(wflow, data = bootstrap_sample)
      
      # Create explainer
      
      metrics <- list(
        logloss = MetricsWeighted::logLoss,
        `ROC AUC` = MetricsWeighted::AUC,
        `% Dev Red` = MetricsWeighted::r_squared_bernoulli
      )
      
      var_names <- names(yhats[[k]]$data)[-1]
      
      pred_fun <- function(m, dat) {
        predict(
          m, dat[, colnames(bootstrap_sample)[-1], drop = FALSE],
          type = "prob"
        )$`.pred_1`
      }
      
      fl <- flashlight(
        model = model_fit,
        label = 'class',
        data = bootstrap_sample,
        y = 'class',
        predict_function = pred_fun,
        metrics = metrics
      )
      
      for (j in seq_along(var_names)) {
        if (ice == TRUE) {
          pd_ <- light_ice(fl, v = paste0(var_names[j]), center = 'first') #ice doesn't work well
        } else {
          pd_ <- light_profile(fl, v = paste0(var_names[j]))
        }
        
        #add number of boostrap.
        bs_rep <- data.frame( bootstrap=rep(i, nrow(pd_$data)))
        
        #response name
        bs_name <- data.frame( response=rep(names(Y[k]), nrow(pd_$data)))
        
        pd_data <-data.frame(cbind(pd_$data),  bs_name, bs_rep) #add bootstrap
        
        pd_raw[[i]][[var_names[j]]] <- pd_data  # Save pd_ as a list element
        
        
      }
    }
    
    return(pd_raw)  # Return pd_raw as a list
  }
  
  bstraps_pd_list <- future_lapply(seq(1, n_response), internal_fit_function, future.seed = TRUE)
  
  return(bstraps_pd_list)
}

# Example usage:
# Replace yhats and other variables with your actual data
# yhats <- list(...)
# result <- mrBootstrap(yhats, num_bootstrap = 10, alpha = 0.05, ice = TRUE)
