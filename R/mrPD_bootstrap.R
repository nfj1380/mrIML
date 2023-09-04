

mrPD_bootstrap <- function(mrBootstrap_obj, vi_obj, X, Y, target, global_top_var = 2) {
  
  n_response <- ncol(Y)
  complete_df <- cbind(Y, X)
  n_data <- ncol(complete_df)
  
  # Internal function to combine objects by name
  bind_rows_by_name <- function(list_obj, object_name) {
    filtered_list <- list_obj[names(list_obj) %in% object_name]
    bind_rows(filtered_list)
  }
  
  internal_fit_function <- function(i) {
    object_name <- names(complete_df[i])
    
    combined_list <- list()  # Create an empty list to store combined objects
    
    for (j in 1:n_response) {
      combined_object <- map_dfr(mrBootstrap_obj[[j]], bind_rows_by_name, object_name)
      
      if (nrow(combined_object) > 0) {
        combined_metadata <- data.frame(target = rep(names(Y)[j], nrow(combined_object)))
        combined_object <- cbind(combined_object, combined_metadata)
        
        combined_list[[j]] <- combined_object   # Append the combined object to the list
      }
    }
    
    combined_df <- do.call(rbind, combined_list)  # Convert the list to a data frame
    
    return(combined_df)
  }
  
  pd_list <- future_lapply(seq_len(n_data), internal_fit_function, future.seed = TRUE)
  
  plot_list <- list()  # Create an empty list to store individual plots
  
  vi_obj <- vi_obj[[1]]  # Extract VI data
  
  G_target_data_avg <- vi_obj %>% 
    filter(response == {{target}}) %>% 
    group_by(var) %>% 
    summarise(mean_imp = mean(sd_value)) %>% 
    arrange(desc(mean_imp))
  
  G_top_vars <- head(G_target_data_avg[order(-G_target_data_avg$mean_imp), ], global_top_var)
  
  # Iterate through each pd_list and create individual plots
  for (k in seq_along(pd_list)) {
    
    df <- pd_list[[k]] %>%
      filter(target == {{target}})
    
    if (names(df)[1] %in% G_top_vars$var) {
      if (is.factor(df[[1]]) || (all(df[[1]] %in% c(0, 1)))) {
        
        d1 <- df %>%
          mutate(class = recode(.[[1]], `0` = "absent", `1` = "present"))
        
        plot <- ggplot(d1, aes(x = class, y = value)) +
          geom_boxplot() +
          labs(x = names(d1)[1], y = paste(target, "prob", sep = " ")) +
          theme_bw()
        
      } else {
        
        d1 <- df %>%
          group_by(bootstrap) %>% 
          rename(class = 1)
        
        plot <- ggplot(d1, aes(x = class, y = value, group = interaction(bootstrap, target)))+
          geom_line(alpha = 0.3) +
          labs(x =  names(df)[1], y = paste(target, "prob", sep = " ")) + 
          theme_bw()
      }
      
      plot_list[[k]] <- plot  # Add the plot to the list
    }
  }
  
  plot_list_updated <- plot_list[sapply(plot_list, function(p) any(p$data$value != 0))]
  
  p <- grid.arrange(grobs = plot_list_updated )
  
  # Create combined plot using the order from G_top_vars
  #combined_plot <- plot_grid(  plot_list_updated =   plot_list_updated[G_top_vars$var], ncol = 1, rel_heights = rep(1, length(G_top_vars$var)))
  
  combined_plot <- plot_grid(p, ncol = 1, rel_heights = rep(1, length(G_top_vars$var)))
  
  return(list(pd_list, combined_plot ))  # Return both pd_list and combined_plot
}
