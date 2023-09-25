#' Calculate covariate partial dependencies for mrIML JSDMs (Joint species distirbution models)
#'
#' This function calculates the covariate partial dependency plot for a specified environmental/host variable.
#' It also filters the taxa based on standard deviation thresholds and visualizes the results.
#'
#' @param yhats A list of model predictions.
#' @param Y The response data.
#' @param X The predictor data.
#' @param X1 Additional predictor data excluding the variable of interest.
#' @param var The variable of interest for calculating the profile.
#' @param sdthresh The standard deviation threshold for filtering taxa (default: 0.05).
#' 
#' @return A plot displaying the covariate profile and change in probability for the specified variable.
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
#'covar <- mr_Covar(yhats, X=X,X1=X1, Y=Y, var='scale.prop.zos', sdthresh =0.01) #sdthrsh just plots taxa responding the most.

#' }

mr_Covar <- function(yhats, Y, X, X1, var,  sdthresh =0.05) {
  
  n_response <- length(Y)
  
  profileData_combined <- lapply(seq_len(n_response), function(i) {
    
    fl <- mrFlashlight(yhats, X = cbind(X, X1[-i]), Y = Y, response = "single", index = i, mode = 'classification')
    var_names <- setdiff(names(yhats[[i]]$data)[-1], names(Y))
    
      pd_data <- light_profile(fl, v = paste(var)) 
      nameCov <- as.data.frame(rep(paste(var), times=nrow(pd_data$data)))
      pd_data_final <- cbind(nameCov, pd_data$data)
      names(pd_data_final)[1] = 'cov'
      names(pd_data_final)[2] = 'cov_grid'
      names(pd_data_final)[4] = paste("value_",i) ###
      
      pd_data_final 
      
  })
  
 
  all_data <- do.call(cbind, profileData_combined)

  plot_data <-  all_data %>% select(all_of(names(all_data)[duplicated(names(all_data))]%>% setdiff('value'))) 
  
  value_data  <- all_data %>% select(contains('value')) 
  colnames(value_data) <- names(Y) 
  
  #filter taxa not strongly responding
  sdev<- value_data %>% 
    summarize(across(everything(), sd)) #double check this works
    
  Xred <- sdev %>%
    select(where(~. >= sdthresh)) %>% t() %>% as.data.frame() %>% 
    rownames_to_column()
  
  names(Xred)[1] <- 'taxa'
  
  r_value_data <- value_data %>%   select(any_of(Xred$taxa))
  
  combi_data_sub <- as.data.frame(cbind(plot_data, r_value_data)) %>% select(-counts, -label, -type)
  
  combi_data_prep_sub <- combi_data_sub %>% 
    pivot_longer(cols = matches(names(Y)), values_to = "value")
  
  p1 <- ggplot(combi_data_prep_sub, aes(x = cov_grid, y = value, colour=name)) +
    geom_line() +
    labs(x = paste0(var), y = "Probability of occurence") +
    scale_color_discrete(name = "Taxa") +
    theme_bw()
  
  #for the complete data
  combi_data<- as.data.frame(cbind(plot_data, value_data)) %>% 
    select(-counts, -label, -type) 
  
  combi_data_prep <-  combi_data %>% 
    pivot_longer(cols = matches(names(Y)), values_to = "value")
  
 
  #cumulation
  # Step 2: Convert "cov_grid" column to factor with specific order if needed. 
  #need to do this for the complete data.
  #change combiprep
  combi_data_prep$cov_grid <- factor(combi_data_prep$cov_grid, levels = unique(combi_data_prep$cov_grid))
  
  # Step 3: Calculate change in values
  df <- combi_data_prep %>%
    arrange(cov_grid) %>%
    group_by(cov_grid) %>%
    mutate(change = ifelse(is.na(value - lag(value)), 0, value - lag(value)))
  
  # Step 4: Plot smoothed change in values
  p2 <-ggplot(df, aes(x = cov_grid, y = change)) +
    geom_boxplot() +
    labs(x = paste0(var), y = "Change in probability") +
    theme_bw()
  
  grid.arrange(p1, p2)

}
