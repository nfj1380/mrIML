#' Bootstrap model predictions
#'
#' This function bootstraps model predictions and generates variable profiles
#' for each response variable based on the provided yhats.
#'
#' @param yhats A list of model predictions (e.g., from mrIMLpredicts).
#' @param num_bootstrap The number of bootstrap samples to generate (default: 10).
#' @param Y The response data.
#' @param mode \code{character}: 'classification' or 'regression' depending on the model type.
#' @param downsample Logical. Should the bootstrap samples be downsampled? (default: FALSE).
#' @return A list containing bootstrap samples of variable profiles for each response variable.
#' @export
#' @examples
#' \dontrun{
#' # Example usage:
#' 
#' # Prepare response data
#' Y <- dplyr::select(Bird.parasites, -scale.prop.zos) %>% 
#'        dplyr::select(sort(names(.))) # Response variables (e.g., SNPs, pathogens, species)
#'
#' # Example list of yhats generated from mrIMLpredicts (assume yhats_rf is defined)
#' yhats_rf <- mrIMLpredicts(...)  # Replace with actual code to generate yhats_rf
#'
#' # Perform bootstrap analysis
#' bs_analysis <- mrBootstrap(yhats = yhats_rf, Y = Y, num_bootstrap = 50, mode = 'classification')
#' }

  
  mrBootstrap <- function(yhats, num_bootstrap = 10, Y=Y, downsample=FALSE, mode='classification') {
    
    n_response <- length(yhats)
    
    pb <- txtProgressBar(min = 0, max = n_response, style = 3) 
    # pb <- progress::progress_bar$new(format = "[:bar] :percent ETA: :eta", total = n_response )
    
    
    internal_fit_function <- function(k) {
      
      setTxtProgressBar(pb, k) #progressbar marker
      # pb$tick() 
      
      features <- colnames(yhats[[k]]$data)[-1]
      
      n <- nrow(yhats[[k]]$data)
      
      pd_raw <- vector("list", num_bootstrap)  # Initialize pd_raw as a list
      
      for (i in 1:num_bootstrap) {
        
        # Initialize the bootstrap sample
        bootstrap_sample <- NULL
        
        
        if (downsample==TRUE) {
          # Determine the number of samples to draw for each class
          class_counts <- table(Y[[k]])
          min_class_count <- min(class_counts)
          sample_size <- min_class_count * length(class_counts)
          
          
          # Sample from each class to balance classes
          for (cls in unique(Y[[k]])) {
            cls_indices <- sample(which(Y[[k]] == cls), size = min_class_count, replace = FALSE)
            bootstrap_sample <- rbind(bootstrap_sample, yhats[[k]]$data[cls_indices, ])
          }
          
        } else {
          
          
          # Convert the data frame to a data table
          data_table <- data.table::as.data.table(yhats[[k]]$data)
          
          # Generate random row indices
          sample_indices <- sample(1:n, replace = TRUE)
          
          # Create the bootstrap sample using data table syntax
          bootstrap_sample <- data_table[sample_indices]
          
        }
        
        # Extract the workflow from the best fit
        wflow <- yhats[[k]]$last_mod_fit %>% tune::extract_workflow()
        
        # Add the bootstrap data to the workflow
        wflow$data <- bootstrap_sample
        
        # Fit the model using the bootstrap sample
        model_fit <- fit(wflow, data = bootstrap_sample)
        
        # Create explainer
        
        var_names <- names(yhats[[k]]$data)[-1]
        
        if(mode=='classification'){
          
          #metric list for flashlight.
          metrics <- list(
            logloss = MetricsWeighted::logLoss,
            `ROC AUC` = MetricsWeighted::AUC,
            `% Dev Red` = MetricsWeighted::r_squared_bernoulli
          )
          
          pred_fun <- function(m, dat) {
            predict(
              m, dat[, colnames(bootstrap_sample)[-1], drop = FALSE],
              type = "prob"
            )$`.pred_1`
          }
        } else {
          
          pred_fun <- function(m, dat) {
            
            predict(m, dat[, colnames(X), drop = FALSE])[[".pred"]]
            
          }
          
          # List of metrics
          
          metrics = list( 
            rmse = MetricsWeighted::rmse,
            `R-squared` = MetricsWeighted::r_squared
          )
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
          
          pd_ <- light_profile(fl, v = paste0(var_names[j]))
          
          #add number of boostrap.
          bs_rep <- data.frame( bootstrap=rep(i, nrow(pd_$data)))
          
          #response name
          bs_name <- data.frame( response=rep(names(Y[k]), nrow(pd_$data)))
          
          pd_data <-data.frame(cbind(pd_$data),  bs_name, bs_rep) #add bootstrap
          
          pd_raw[[i]][[var_names[j]]] <- pd_data  # Save pd_ as a list element
          
        }
      }
      
      gc() #clear junk
      
      return(pd_raw)  # Return pd_raw as a list
    }
    
    bstraps_pd_list <- future_lapply(seq(1, n_response), internal_fit_function, future.seed = TRUE)
   
     return(bstraps_pd_list) 
}
    
  
