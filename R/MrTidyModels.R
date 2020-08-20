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

mrIMLpredicts<- function(X, Y, model1, balance_data ='no') { 
  
  n_response<- length(X)
  # Run model 1 for each parasite; a simple logistic regression with a single covariate
  # in this case but note that model 1 can be any model of the user's choice, 
  # from simple regressions to complex hierarchical or deep learning models.
  # Different structures can also be used for each species to handle mixed multivariate outcomes
  
  #options(show.error.messages= FALSE) not working
  #sink(type="message")
  
  mod1_perf <- NULL #place to save performance matrix
  
  #yhats <- for(i in 1:length(X)) {
  yhats <- lapply(seq(1,n_response), function(i){
    #rhats <- lapply(seq(1, n_variables), function(species){
    
    #not needed for this model
    #OtherSNPs <- as.data.frame(X[-1]) 
    #OtherSNPs[OtherSNPs=='Negative'] <- 0
    #OtherSNPs[OtherSNPs== 'Positive'] <- 1 #could do a PCA/PCoA?
    #OtherSNPsa <-apply(OtherSNPs, 2, as.numeric) 
    
    data <- cbind(X[i ], Y) 
    colnames(data)[1] <- c('class') #define response variable
    
    data$class<- as.factor(data$class)
    
    #data<-data[complete.cases(data)] #removes NAs but there must be a conflict somewhere
    data_split <- initial_split(data, prop = 0.75)
    
    # extract training and testing sets
    data_train <- training(data_split)
    data_test <- testing(data_split)
    
    #10 fold cross validation
    data_cv <- vfold_cv(data_train, v= 10) 
      
    if(balance_data == 'down'){ 
      data_recipe <- training(data_split) %>%
        recipe(class ~., data= data_train) %>% #if downsampling is needed
        themis::step_downsample(class)
      }
    
    if(balance_data == 'up'){ 
      data_recipe <- training(data_split) %>%
        recipe(class ~., data= data_train) %>% #if downsampling is needed
        themis::step_upsample(class)
    }
    
    if(balance_data == 'rose'){
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
    
    # Fit model one for each parasite; can easily modify this so that the user
    # can specify the formula necessary for each species as a list of formulas
    
    # Find best tuned model
    k3_20 <- expand_grid(neighbors = 3:20)
    res_tune <-
      mod_workflow %>%
      tune::tune_grid(resamples = data_cv,
                      grid = k3_20,
                      metrics = yardstick::metric_set(roc_auc),#remove accuracy
                      control= tune:: control_resamples(save_pred = TRUE),
                      collect_metrics()) #this isnt working properly - what is neighbors in this context?
    
    # get the results of tunning
    tune_results <- res_tune %>% 
      show_best(metric = "roc_auc", n = 5)
    
    param_final <- res_tune %>%
      select_best(metric = "roc_auc")
    
    #update model
    final_workflow <- mod_workflow %>%
      finalize_workflow(param_final)
    
    mod1_k <- final_workflow %>%
      fit(data = data_train)
    
    # the best model fit
   # set.seed(345)
    best_mod_fit <- 
     mod_workflow %>% 
     last_fit(data_split,res_tune)
      select_best(res_tune, metric = "roc_auc")
    
   # best_mod_fit$.workflow
  
## add the best for the tuning  option for user (Gustavo)# 
## turn it off if you want 
    
    
    #fit on the training set and evaluate on test set. Not needed 
    #last_fit(data_split) 
    
    # Calculate probability predictions for the fitted training data. 
    
    yhatO <- predict(mod1_k, new_data = data_train, type='prob' )
    
    yhat <- yhatO$.pred_1
    
    #' Calculate deviance residuals 
    resid <- devianceResids(yhatO, data_train$class)
    
    #predictions based on testing data
    yhatT <- predict(mod1_k, new_data = data_test, type='prob') %>% 
      bind_cols(data_test %>%
                  select(class))
  
    })
    
    
    list(mod1_k = mod1_k, res_tune=res_tune, best_mod_fit=best_mod_fit, tune_results=tune_results, data=data, data_testa=data_test, data_train=data_train, yhat = yhat, yhatT = yhatT, resid = resid) 
  } 


