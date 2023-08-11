mrStacksVI <- function(yhats_stacks, Y = Y, ice = FALSE, local_top_var=5) {
  
  n_response <- length(yhats_stacks)
  
  pb <- txtProgressBar(min = 0, max = n_response, style = 3) 
  
  internal_fit_function <- function(k) {
    
    setTxtProgressBar(pb, k) #progressbar marker
    
    features <- colnames(yhats_stacks[[k]]$data)[-1]
    
    n <- nrow(yhats_stacks[[k]]$data)
    
    # Extract the workflow from the best fit
    wflow <- yhats_stacks[[k]]$last_mod_fit %>% extract_workflow()
    
    # Fit the model using the bootstrap sample
    model_fit <- fit(wflow, data = yhats_stacks[[k]]$data)
    
    # Prediction function for regression
    pred_fun <- function(m, dat) {
      predict(m, dat[, names(dat)[-1], drop = FALSE])[[".pred"]]
    }
    
    # List of metrics
    metrics <- list( 
      rmse = MetricsWeighted::rmse,
      `R-squared` = MetricsWeighted::r_squared
    )
    
    var_names <- names(yhats_stacks[[k]]$data)[-1]
    
    fl <- flashlight(
      model = model_fit,
      label = 'class',
      data = yhats_stacks[[k]]$data,
      y = 'class',
      predict_function = pred_fun,
      metrics = metrics
    )
    
    pd_raw <- vector("list", length = length(var_names))  # Initialize pd_raw list
    
    for (j in seq_along(var_names)) {
      
      if (ice == TRUE) {
        pd_ <- light_ice(fl, v = paste0(var_names[j]), center = 'first') # ice doesn't work well
      } else {
        pd_ <- light_profile(fl, v = paste0(var_names[j]))
      }
      # Add taxa
      taxa_names <- data.frame(response = rep(names(Y[k]), nrow(pd_$data)))
      feature_names <- data.frame(feature = rep(names(pd_$data)[1], nrow(pd_$data)))
      
      names(pd_$data)[1] <- 'feature' #make sure a merge can happen
      pd_data <- data.frame(cbind(feature_names, pd_$data, taxa_names))
      
      pd_raw[[j]] <- pd_data  # Store data frame in pd_raw list
    }
    
    return(pd_raw)  # Return the pd_raw list for this response variable
  }
  
  pd_list <- lapply(seq_along(yhats_stacks), internal_fit_function)
  
  #put it all together
  merged_data <- bind_rows(pd_list)
  
  result <- merged_data %>%
    group_by(response, feature) %>%
    summarise(sd_value = sd(value))%>% 
    ungroup() 
  
  p1 <- ggplot(result, aes(y=reorder(feature,sd_value), x=sd_value))+
    geom_boxplot()+
    theme_bw()+
    labs(x="Importance", y="Features")+
    theme(
      axis.title.x = element_text(size = 8),  # Adjust the size as needed
      axis.title.y = element_text(size = 8))
  
  filtered_data <-   result %>%
    group_by(response) %>%
    summarise(sd_value = mean(sd_value)) %>%
    right_join(ModelPerf[[1]], by = join_by(response)) %>% 
    filter( rsquared > threshold)
  
  #plot_boxplots <- function(vi_df){  
  # Create a list to store the plots
  plot_list <- list()
  
  # Generate boxplots for each target
  for (l in 1:nrow(filtered_data)) {
    
    target <- filtered_data$response[l]###
    
    # Filter the data for the current target
    target_data <-  merged_data %>%
      filter(response == {{target}})
    
    target_data_avg <- target_data %>% 
      group_by(feature) %>% 
      summarise(mean_imp = sd(value)) 
    
    # Create the boxplot for the target
    
    #get top variables
    top_vars <- head(target_data_avg[order(-target_data_avg$mean_imp), ], local_top_var)
    
    target_data_final_df <- target_data %>% 
      filter(feature %in% top_vars$feature) %>% 
      group_by(feature) %>% 
      summarise(sd_value = sd(value)) 
      
    plot <- ggplot(target_data_final_df, aes(x = sd_value, y = reorder(feature , sd_value))) +
      geom_bar(stat = "identity") +
      theme_bw()+
      labs(x="Importance", y=target)+
      scale_y_discrete(label=abbreviate)+
      theme(
        axis.title.x = element_text(size = 8),  # Adjust the size as needed
        axis.title.y = element_text(size = 8))+
      coord_flip() 
    
    # Add the plot to the list
    plot_list[[l]] <- plot
  }
  
  
  # Arrange and display the plots in a grid
  p2 <- grid.arrange(grobs = plot_list)
  
  #publication ready plot  
  combined_plot <- plot_grid(p1, p2, rel_heights = c(1, 0.5), labels = "auto")
  
  return(list(merged_data,  combined_plot))
  
}
