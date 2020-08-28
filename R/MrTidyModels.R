#'mrIMLpredicts:  Wrapper to generate multi-response tidy models.This function produces yhats that are used in the stack function as well 
#'as saving all model characteristics for subsequent functions.
#'
#'This function fits separate classication models for each response variable in a dataset. 
#'@param Y A \code{dataframe} is a response variable data set (species, OTUs, SNPs etc).
#'@param X A \code{dataframe} represents predictor or feature data.
#'@param balance_data A \code{character} 'up', 'down' or 'no'. 

#'@param Model 1 A \code{list} can be any model from the tidy model package. See examples.
#'
#'@examples
#'model1 <- #model used to generate yhat
#'specify that the model is a random forest
#'logistic_reg() %>%
#' # select the engine/package that underlies the model
#'set_engine("glm") %>%
#'  # choose either the continuous regression or binary classification mode
#'  set_mode("classification")
#'  
#'@details Y (response variables) should be binary (0/1). Rows in X (features) have the same id (host/site/population)
#'  as Y. 
#'  Class imblanace can be a real issue for classification analyses. Class imbalance can be addressed for each
#'   response variable using 'up' (upsampling using ROSE bootstrapping), 'down' (downsampling) 
#'or 'no' (no balancing of classes).


mrIMLpredicts<- function(X, Y, model1, balance_data ='up') { 
  
  n_response<- length(X)
  # Run model 1 for each parasite; a simple logistic regression with a single covariate
  # in this case but note that model 1 can be any model of the user's choice, 
  # from simple regressions to complex hierarchical or deep learning models.
  # Different structures can also be used for each species to handle mixed multivariate outcomes
  
  mod1_perf <- NULL #place to save performance matrix
  
  #yhats <- for(i in 1:length(X)) {
  yhats <- lapply(seq(1,n_response), function(i){
    #rhats <- lapply(seq(1, n_variables), function(species){
    
    #not needed for this model
    #OtherSNPs <- as.data.frame(X[-1]) 
    #OtherSNPs[OtherSNPs=='Negative'] <- 0
    #OtherSNPs[OtherSNPs== 'Positive'] <- 1 #could do a PCA/PCoA?
    #OtherSNPsa <-apply(OtherSNPs, 2, as.numeric) 
    
    data <- cbind(X[i], Y) 
    colnames(data)[1] <- c('class') #define response variable
    
    data$class<- as.factor(data$class)
    
    #data<-data[complete.cases(data)] #removes NAs but there must be a conflict somewhere
    set.seed(100)
    data_split <- initial_split(data, prop = 0.75)
    #data_splitalt <- initial_split(data, strata = class)
    
    # extract training and testing sets
    data_train <- training(data_split)
    data_test <- testing(data_split)
    
    # extract training and testing sets stata
    #data_trainalt <- training(data_splitalt)
    #data_testalt <- testing(data_splitalt)
    
    #10 fold cross validation
    data_cv <- vfold_cv(data_train, v= 10) #not used yet.
    
    if(balance_data == 'down'){ 
      data_recipe <- training(data_split) %>%
        recipe(class ~., data= data_train) %>% #if downsampling is needed
        themis::step_downsample(class)
    }
    
    if(balance_data == 'up'){
      data_recipe <- training(data_split) %>%
        recipe(class ~., data= data_train) %>%
        themis::step_rose(class) #ROSE works better on smaller data sets. SMOTE is an option too.
    }
    
    
    if(balance_data == 'no'){ 
      data_recipe <- training(data_split) %>% #data imbalance not corrected 
        recipe(class ~., data= data_train)
    }
    
    #optional recipe ingredients
    #step_corr(all_predictors()) %>% # removes all corrleated features
    #step_center(all_predictors(), -all_outcomes()) %>% #center features
    #step_scale(all_predictors(), -all_outcomes()) %>% #scale features
    
   
    
    mod_workflow <- workflow() %>%
      # add the recipe
      add_recipe(data_recipe) %>%
      # add the model
      add_model(model1)

    
     ## full tunning 
    
     tune_m<-tune::tune_grid(mod_workflow,
                       resamples = data_cv,
                            grid = 10)
     
     # select the best model
      best_m <- tune_m %>%
       select_best("roc_auc")
      
      # final workflow
      final_model <- 
        mod_workflow %>% 
       finalize_workflow(best_m)
      
    # Fit model one for each parasite; can easily modify this so that the user
    # can specify the formula necessary for each species as a list of formulas
     mod1_k <- final_model %>%
      fit(data = data_train)
    
    #mod1_k %>%
      #fit_resamples(resamples = data_cv)

# keep the tune all list 
    
    # the last fit
    set.seed(345)
    last_mod_fit <- 
      final_model %>% 
      last_fit(data_split)
    
    #fit on the training set and evaluate on test set. Not needed 
    #last_fit(data_split) 
    
    # Calculate probability predictions for the fitted training data. 
    
    yhatO <- predict(mod1_k, new_data = data_train, type='prob' )
    
    yhat <- yhatO$.pred_1
    
    #' Calculate deviance residuals 
    resid <- devianceResids(yhatO, data_train$class )
    
    #predictions based on testing data
    yhatT <- predict(mod1_k, new_data = data_test, type='class' ) %>% 
      bind_cols(data_test %>% select(class))
    
    
   list(mod1_k = mod1_k, last_mod_fit=last_mod_fit,tune_m=tune_m, data=data, data_testa=data_test, data_train=data_train, yhat = yhat, yhatT = yhatT, resid = resid) 
  })
}
