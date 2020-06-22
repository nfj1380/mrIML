#'mrTidyModels:  Wrapper to generate multi-response tidy models. 
#'
#'
#'This function fits separate classication models for each response variable in a dataset. 
#'@param X  #is a response variable data set (specie, SNPs etc).
#'@param Y #represents predictor or feature data.
#'Row id's must match bewtween X and Y. Will add something to make sure of this.

#'@param Model 1 #can be any model from the tidy model package in the follwoing format
#'model1 <- 

#logistic_reg() %>%
#'set_engine("glm") %>%
#'set_mode("classification")
#'
#'
#'This function produces yhats that are used in the stack function as well 
#'as saving all model characteristics for subsequent functions.

mrTidyPredict<- function(X, Y, model1) {
  
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
    colnames(data)[1] <- c('class')
    
    data$class<- as.factor(data$class)
    
    #data<-data[complete.cases(data)] #removes NAs but there must be a conflict somewhere
    data_split <- initial_split(data, prop = 0.75)
    
    # extract training and testing sets
    data_train <- training(data_split)
    data_test <- testing(data_split)
    #10 fold cross validation
    data_cv <- vfold_cv(data_train, v= 10)
    
    data_recipe <- training(data_split) %>%
      recipe(class ~., data= data_train)
    
    mod_workflow <- workflow() %>%
      # add the recipe
      add_recipe(data_recipe) %>%
      # add the model
      add_model(model1)
    
    # Fit model one for each parasite; can easily modify this so that the user
    # can specify the formula necessary for each species as a list of formulas
    
    mod1_k <- mod_workflow %>%
      fit(data = data_train)
    
    #fit on the training set and evaluate on test set. Not needed 
    #last_fit(data_split) 
    
    # Calculate probability predictions for the fitted training data. The steps belwo didnt work
    
    yhatO <- predict(mod1_k, new_data = data_train, type='prob' )
    
    yhat <- yhatO$.pred_1
    
    #' Calculate deviance residuals 
    resid <- devianceResids(yhatO, data_train$class )
    
    #predictions based on testing data
    yhatT <- predict(mod1_k, new_data = data_test, type='class' ) %>% 
      bind_cols(data_test %>% select(class))
    
    
    list(mod1_k = mod1_k, yhat = yhat,yhatT = yhatT, resid = resid) 
  })  

  
}

