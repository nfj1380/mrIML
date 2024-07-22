#' Plot Model Performance Comparison
#'
#' Create visualizations to compare the performance of two models based on their performance metrics.
#'
#' @param ModelPerf1 Dataframe of model performance metrics for the first model to compare.
#' @param ModelPerf2 Dataframe of model performance metrics for the second model to compare.
#' @param mod_names Character vector of model names. Default is \code{c('combined', 'Xonly_model')}.
#' @param mode Character string indicating whether the mode is 'classification' or 'regression'. Default is 'classification'.
#' @return A list containing:
#' \item{p1}{A ggplot object for the boxplot of model performance metrics.}
#' \item{p2}{A ggplot object for the barplot of differences in performance metrics.}
#' \item{wide_df}{A dataframe with the wide format of model performance metrics and their differences.}
#' @export
#' @examples
#'plots <- mrPerformancePlot(ModelPerf1 =ModelPerf_lm, ModelPerf2 = ModelPerf_rf, mod_names=c('linear_reg','rand_forest'), mode='regression' )
#' 
mrPerformancePlot <- function(ModelPerf1 = NULL, ModelPerf2 = NULL,
                              mod_names=c('combined','Xonly_model'), mode='classification' ) {
  
  m1 <- ModelPerf1[[1]]
  m2 <- ModelPerf2[[1]]
  
  d_length <- nrow(m1) #has to be the same
  
  if (mode=='classification'){
    
    m1$mcc <- as.numeric(m1$mcc) #extract mcc
    
    m1 <- m1 %>%
      dplyr::rename(metric = mcc)
    
    model1_df <- m1 %>% tibble(model_type = rep(mod_names[1], d_length)) %>% 
      replace_na(list(metric = 0)) #NAs happen when there is no variance in the model predictions (i.e.all 0s), for these make NAs=0
    
    m2$mcc <- as.numeric(m2$mcc)
    
    model2_df <- m2 %>% tibble(model_type = rep(mod_names[2], d_length)) %>% 
      dplyr::rename(metric=mcc) %>%
      replace_na(list(metric = 0))
 
  } else {
    
    model1_df <- m1 %>% tibble(model_type = rep(mod_names[1], d_length)) %>% 
      dplyr::rename(metric=rmse)
    
    model2_df <- m2 %>% tibble(model_type = rep(mod_names[2], d_length)) %>% 
      dplyr::rename(metric=rmse)
    
  }
    model_compare_df <- bind_rows(model1_df , model2_df )  
  
#detect outliers function
findoutlier <- function(x) {
  return(x < quantile(x, .25) - 1.5*IQR(x) | x > quantile(x, .75) + 1.5*IQR(x))
}  
  #create data frame
  plot_df <- model_compare_df %>%
    group_by(model_name) %>%
    mutate(outlier = ifelse(findoutlier(metric), metric, NA))
   

  if(mode=='classification') {
    
  p1 <- ggplot(plot_df, aes(x=model_name, y=metric)) +
    geom_boxplot() +
    geom_text(aes(label=outlier), na.rm=TRUE, hjust=-.5) +
    theme_bw()+
    labs(y='MCC')

  
  } else {
    
  p1 <- ggplot(plot_df, aes(x=model_name, y=metric)) +
    geom_boxplot() +
    geom_text(aes(label=outlier), na.rm=TRUE, hjust=-.5) +
    theme_bw()+
    labs(y='RMSE')
  }
  
  #individual taxa 
  
  wide_df <- plot_df %>%
    select(response, model_name, metric, outlier) %>% 
    pivot_wider(names_from = model_name, values_from = metric)
  
  # Calculate differences from the first model (combined)
  
  wide_df <- wide_df %>%
    mutate( diff_mod1_2 = .[[4]] - .[[3]])
  
  
  # Reshape back to long format for plotting
  long_df <- wide_df %>%
    pivot_longer(cols = starts_with("diff_"), names_to = "comparison", values_to = "difference")
  
  # Plot the differences
  
  if(mode=='classification'){
    
  p2 <- ggplot(long_df, aes(x=reorder(response, difference, decreasing =T), y=  difference))+
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(x = "Response", y = "Difference in MCC", title = paste(mod_names[2]," vs ",paste(mod_names[1])) ) +
    theme_bw()
  
  } else {
    
    p2 <- ggplot(long_df, aes(x=reorder(response, difference, decreasing =T), y=  difference))+
      geom_bar(stat = "identity") +
      coord_flip() +
      labs(x = "Response", y = "Difference in RMSE", title = paste(mod_names[2]," vs ",paste(mod_names[1])) ) +
      theme_bw()
    
  }

  return(list(p1, p2, wide_df))
}
