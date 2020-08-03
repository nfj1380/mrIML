#'Generates covariance matrix
#'
#'Problem is inf values in yhats using deviance function

response_covariance = function(yhats, covariance_mod){
  
  #test data is now within the yhat object
  #test.data <- yhats[[1]]$data_testa
  
  # Determine the total number of outcome variables in the yhat data
  n_variables <- length(yhats)
  
  # For each outcome, use model 1 residuals (resid from the yhats list) as outcomes and 
  # predictions from each other outcome's mod 1 (yhat) as predictors to learn the covariance matrix
  cat('Fitting sequential', covariance_mod$engine, 'models to learn the multivariate covariance matrix...\n')
  rhats <- lapply(seq(1, n_variables), function(response){
    outcome <- yhats[[response]]$resid #response 
    predictors <- do.call(cbind, purrr::map(yhats, 'yhat'))[, -response]#response
    colnames(predictors) <- paste0('pred_', seq(1, ncol(predictors)))
    
    data_R <-  cbind(data.frame(Y = outcome),
                     as.matrix(predictors))
    
    data_split_R <- initial_split(data_R, prop = 0.75) #75% training and 25 %testing data
    
    data_train_R <- training(data_split_R) 
    data_test_R <- testing(data_split_R)
    
    # extract training and testing sets
    
    #10 fold cross validation
    #data_cv <- vfold_cv(data_train_R, v= 10) #not properly implemented yet
    
    formula <- as.formula(paste0("Y~", 
                                 paste0(colnames(predictors),collapse = "+")))
    
    #data_recipe_R <- training(data_R) %>% #data imbalance not corrected 
     #
    data_recipe_R <- recipe(formula, data= data_R)
    
    
    mod2_workflow <- workflow() %>%
      # add the recipe
      add_recipe(data_recipe_R) %>%
      # add the model
      add_model(covariance_mod)
    
    set.seed(123)
    #10 fold CV to check performance
    folds <- vfold_cv( data_R,  v = 10) # the complete data?
    flights_fit_rs <- fit_resamples(mod2_workflow, folds)
    collect_metrics(flights_fit_rs)# need to incorporate this.
    
    # Fit model one for each parasite; can easily modify this so that the user
    # can specify the formula necessary for each species as a list of formulas
    
    mod2 <- mod2_workflow %>%
      fit(data =  data_R )
    
    # the last fit
    #set.seed(345)
    last_mod2_fit <- 
      mod2_workflow %>% 
      last_fit(data_split_R) %>% 
      purrr::pluck(".workflow",1) %>%   
      pull_workflow_fit()
    
    
    
    # Return the model as well as the 'multivariate residual adjustments'
    list(rhat_k = predict(mod2, new_data = select(data_R, -Y)), mod2 = last_mod2_fit , pred_names = colnames(predictors))
  })
}
  