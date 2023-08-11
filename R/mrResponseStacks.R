


mrResponseStacks<- function (yhats, Model, taxa=Y, prop=0.5, tune_grid_size=5, k=10, seed=1234){
  
  #gather residuals
  dList <- yhats %>% purrr::map(pluck("deviance"))
  deviance_df <- data.frame(do.call(cbind, dList))
  deviance_df_noinf<- deviance_df %>% mutate_all(~ ifelse(. == -Inf, -3, .))
  deviance_df_noinf<- deviance_df_noinf %>% mutate_all(~ ifelse(. == Inf, 3, .))
  # glimpse(deviance_df_noinf)
  
  final_deviance <- deviance_df_noinf
  
  n_response<- length(Y)
  
  mod1_perf <- NULL #place to save performance matrix
  
  internal_fit_function <- function( i ){
    
    data <- cbind(final_deviance[i], taxa[,-i]) ###
    colnames(data)[1] <- c('class') #define response variable 
    
    set.seed(seed)
    
    data_split <- initial_split(data, prop = prop)
    
    # extract training and testing sets
    data_train <- training(data_split)
    data_test <- testing(data_split)
    
    #n fold cross validation
    data_cv <- vfold_cv(data_train, v= k) 
    
    data_recipe <- training(data_split) %>% 
      recipe(class ~., data= data_train)
    
    
    mod_workflow <- workflow() %>%
      # add the recipe
      add_recipe(data_recipe) %>%
      # add the model
      add_model(Model)
    
    ## Tune the model
    
    tune_m<-tune::tune_grid(mod_workflow,
                            resamples = data_cv,
                            grid = tune_grid_size) 
    
    
    # select the best model
    best_m <- tune_m %>%
      select_best("rmse")
    
    
    # final model specification
    final_model <- finalize_workflow(mod_workflow,
                                     best_m )
    # now to fit the model
    mod1_k <- final_model %>%
      fit(data = data_train)
    
    yhatO <- predict(mod1_k, new_data = data_train ) 
    
    yhat <- yhatO$.pred
    
    #predictions based on testing data
    yhatT <- predict(mod1_k, new_data = data_test) # %>% 
    #bind_cols(data_test %>% select(class))
    
    deviance=NULL
    
    # the last fit. Useful for some functionality
    
    last_mod_fit <- 
      final_model %>% 
      last_fit(data_split)
    
    #save data
    list(mod1_k = mod1_k, last_mod_fit=last_mod_fit,tune_m=tune_m, data=data, data_testa=data_test, data_train=data_train, yhat = yhat)
    
    
  }
  
  yhats <- future_lapply(seq(1,n_response), internal_fit_function, future.seed = TRUE)
  
}