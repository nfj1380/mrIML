#' Calculates and plots bootstrapped variable importance
#' @param mrBootstrap_obj The object containing model bootstrapping results.
#' @param ModelPerf A list containing model performance metrics for each response variable.
#' @param X The predictor data.
#' @param Y The response data.
#' @param threshold The threshold for model performance (AUC) below which variables are filtered (default: 0.1).
#' @param global_top_var The number of top global variables to display (default: 10).
#' @param local_top_var The number of top local variables for each response to display (default: 5).
#' 
#' @return A list containing variable importance data and a combined plot.
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' #set up analysis
#' Y <- dplyr::select(Bird.parasites, -scale.prop.zos)%>% 
#' dplyr::select(sort(names(.)))#response variables eg. SNPs, pathogens, species....
#' X <- dplyr::select(Bird.parasites, scale.prop.zos) # feature set

#' X1 <- Y %>%
#' dplyr::select(sort(names(.)))
#'model_rf <- 
#' rand_forest(trees = 100, mode = "classification", mtry = tune(), min_n = tune()) %>% #100 trees are set for brevity. Aim to start with 1000
#' set_engine("randomForest")
#' yhats_rf <- mrIMLpredicts(X=X, Y=Y,
#'X1=X1,'Model=model_rf , 
#'balance_data='no',mode='classification',
#'tune_grid_size=5,seed = sample.int(1e8, 1),'morans=F,
#'prop=0.7, k=5, racing=T) #
#'ModelPerf <- mrIMLperformance(yhats_rf, Model=model_rf, Y=Y, mode='classification')
#'
#'bs_analysis <- mrBootstrap(yhats=yhats_rf,Y=Y, num_bootstrap = 5)
#'bs_impVI <- mrVI_bootstrap(mrBootstrap_obj=bs_analysis, ModelPerf=ModelPerf, 
#'threshold=0.8,  X=X, Y=Y, global_top_var=10,
#'local_top_var=5)
#'} 
#' 
mrVI_bootstrap <- function(mrBootstrap_obj, ModelPerf, X, Y,
                           threshold=0.1, global_top_var=10, local_top_var=5) {
  #threshold is AUC now (less NAs)
  #this stops an annoying message from showing up
  options(dplyr.summarise.inform = FALSE)
  
  n_response <- ncol(Y)
  complete_df <- cbind(Y, X)
  n_data <- ncol(complete_df)
  
  internal_fit_function <- function(i) {
    object_name <- names(complete_df[i])
    
    bind_rows_by_name <- function(list_obj, object_name) {
      filtered_list <- list_obj[names(list_obj) %in% object_name]
      bind_rows(filtered_list)
    }
    
    combined_list <- list()  # Create an empty list to store combined objects
    
    for (j in 1:n_response) {
      combined_object <- map_dfr(mrBootstrap_obj[[j]], bind_rows_by_name, object_name)
      
      if (nrow(combined_object) > 0) {
        
        combined_objfinal <- combined_object #this makes sure it doesnt save the duplicates
        
        combined_list[[j]] <-  combined_objfinal  # Append the combined object to the list
      }
    }
    
    combined_df <- do.call(rbind, combined_list)  # Convert the list to a data frame
    
    return(combined_df)
  }
  
  pd_list <- future_lapply(seq_len(n_data), internal_fit_function, future.seed = TRUE)
  
  
  vi_list <- list()  # Create an empty list to store the plots
  
  for (i in seq_along(pd_list)) {
    
    #extract data for each variable in the models
    df <-  pd_list[[i]]
    
    #calculate the standard deviation of each bootstrap for each taxa
    result <- df %>%
      group_by(response, bootstrap) %>%
      summarise(sd_value = sd(value))%>% 
      ungroup() 
    
    #add the variables names back in
    
    df_names <- rep(names(df[1]), nrow(result))
    
    result_update <- cbind(result, var=df_names)
    
    vi_list[[i]] <- result_update
    
  }
  
  vi_df <- do.call(rbind, vi_list)
  #NB#could have a groupCOv argument here
  
  G_target_data_avg <- vi_df %>% 
    group_by(var) %>% 
    summarise(mean_imp = mean(sd_value)) 
  
  # Create the boxplot for the target
  
  #get top variables
  G_top_vars <- head(G_target_data_avg[order(-G_target_data_avg$mean_imp), ], global_top_var)
  
  G_target_data_final_df <- vi_df %>% 
    filter(var %in% G_top_vars$var)
  
  #plot overall importance 
  p1 <- ggplot(G_target_data_final_df, aes(y=reorder(var,sd_value), x=sd_value))+
    geom_boxplot()+
    theme_bw()+
    labs(x="Importance", y="Features")+
    theme(
      axis.title.x = element_text(size = 8),  # Adjust the size as needed
      axis.title.y = element_text(size = 8)
    )
  
  grid.arrange(p1)
  
  #boxplots for each individual predictor
  
  # Filter the data based on the threshold and calculate mean values for each target
  
  filtered_data <-  vi_df %>%
    group_by(response) %>%
    summarise(sd_value = mean(sd_value)) %>%
    right_join(ModelPerf[[1]], by = join_by(response)) %>% 
    filter(roc_AUC > threshold) %>% 
    ungroup()
  
  #plot_boxplots <- function(vi_df){  
  # Create a list to store the plots
  plot_list <- list()
  
  # Generate boxplots for each target
  for (k in 1:nrow(filtered_data)) {
    
    target <- filtered_data$response[k]###
    
    # Filter the data for the current target
    target_data <-  vi_df %>%
      filter(response == {{target}})
    
    target_data_avg <- target_data %>% 
      group_by(var) %>% 
      summarise(mean_imp = mean(sd_value)) %>% 
      ungroup()
    
    # Create the boxplot for the target
    
    #get top variables
    top_vars <- head(target_data_avg[order(-target_data_avg$mean_imp), ], local_top_var)
    
    target_data_final_df <- target_data %>% 
      filter(var %in% top_vars$var)
    
    plot <- ggplot(target_data_final_df, aes(x = sd_value, y = reorder(var, sd_value))) +
      geom_boxplot() +
      theme_bw()+
      labs(x="Importance", y=target)+
      scale_y_discrete(label=abbreviate)+
      theme(
        axis.title.x = element_text(size = 8),  # Adjust the size as needed
        axis.title.y = element_text(size = 8)
      )
    
    # Add the plot to the list
    plot_list[[k]] <- plot
  }
  
  
  # Arrange and display the plots in a grid
  p2 <- grid.arrange(grobs = plot_list)
  
  #publication ready plot  
  combined_plot <- plot_grid(p1, p2, rel_heights = c(1, 0.5), labels = "auto")
  
  return(list(vi_df, combined_plot))
  
}

