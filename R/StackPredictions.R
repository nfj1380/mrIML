#'Generate predictions for multiple outcomes using stacked multivariate learning
#'
#'
#'This function fits separate penalized models to the residuals for each outcome variable 
#'to estimate parameters of complex multivariate variance covariance matrices. 
#'The initial predictions (supplied by the user) are then updated using the learned
#'covariance parameters.
#'
#'@param yhats A \code{list} with one slot per binary outcome variable. Each outcome's slot
#'should be a \code{list} containing the following named items:
#'\itemize{
#'    \item \code{mod_1}: The full model object used to model the binary outcome variable
#'    \item \code{yhat}: Predictions from \code{mod_1} on the outcome (probability) scale
#'    \item \code{resid}: Residuals from \code{mod_1}. For binary outcomes,
#'    these residuals should be deviance residuals (see examples for how to calculate these)
#'    }
#'@param data.test A \code{dataframe} containing the testing split for cross-validation. This object
#'should have columns for all covariates that were used in the \code{mod_1}'s for generating 
#'out-of-sample predictions
#'@param covariance_mod A \code{character} string specifying which kind of algorithm to use for
#'learning the covariance matrix. Allowed options are currently 'gbm', 'gam' or 'lm'. Note that
#'for gbm models, \code{1000} regression trees will be stored for each outcome variable. Due to memory constraints, 
#'datasets with \code{>50} outcomes therefore cannot use gbm and will instead use a penalized gam by default
#'@return A \code{list} of length \code{length(yhats)} containing:
#'\itemize{
#'    \item \code{probability_preds}: Predictions for \code{data.test} on the probability scale
#'    \item \code{binary_preds}: Predictions for \code{data.test} on the binary scale, where
#'    \code{probability_preds >= 0.5 <- 1}
#'    }
#'@references Xing, L, Lesperance, ML and Zhang, X (2020). Simultaneous prediction of multiple outcomes 
#'using revised stacking algorithms. Bioinformatics, 36, 65-72.
#'
#'@details Separate algorithms are fit to the residuals of each outcome, using the predictions
#'of the remaining outcomes (taken from the \code{yhat} slots in the list \code{yhats}) as covariates.
#'The ability to use complex interactions, such as can be learned in a \code{gbm} model, allows 
#'covariance matrices to be approximated and used to update existing out-of-sample predictions for 
#'each outcome. Crucially, this can be done without already knowing the out-of-sample 'truths' for 
#'each outcome variable, overcoming a major limitation of common joint species distribution models (where the other
#'species' out-of-sample occurrences must be known before generating predictions for the focal species)
#'
#'@examples
#'if(!require(MRFcov)){
#'devtools::install_github('nicholasjclark/MRFcov')}
#'
#'# Load the MRFcov library to access the testing dataset
#'library(MRFcov)

#could merge this into one function i.e. model 1 and model 2 combined?

data("Bird.parasites")
#'
#'
if(!require(devtools)){
install.packages('devtools')
}

if(!require(caret)){
  install.packages('caret')
}

if(!require(gbm)){
  install.packages('gbm')
}

if(!require(MRFcov)){
  devtools::install_github('nicholasjclark/MRFcov')
}

if(!require(tidymodels)){
  install.packages('tidymodels')
}


# Load the MRFcov library to access the testing dataset
library(MRFcov)
#'# A simple split of the data into training and testing subsets
#'train.dat <- Bird.parasites[1:350, ]
#'test.dat <- Bird.parasites[351:449, ]
#'
#'# Run model 1 for each parasite; a simple logistic regression with a single covariate
#'# in this case but note that model 1 can be any model of the user's choice, 
#'# from simple regressions to complex hierarchical or deep learning models.
#'
#'# Different structures can also be used for each species to handle mixed multivariate outcomes


X <- select(Bird.parasites, -scale.prop.zos) #response variables eg. SNPs, pathogens, species....
Y <- select(Bird.parasites, scale.prop.zos) # feature set



#user can specify any tidy model here. 

#user provides the models they would like for each component (model1, model2 for the second bit)
model1 <- #model used to generate yhat
  # specify that the model is a random forest
  logistic_reg() %>%
  # select the engine/package that underlies the model
  set_engine("glm") %>%
  # choose either the continuous regression or binary classification mode
  set_mode("classification")

#----------------------------------------------
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
    
    #calculate metrics on testing data
    mathews <- mcc(yhatT,class, .pred_class)
    mathews <- mathews$.estimate
    sen <- sens(yhatT,class, .pred_class)
    sen <- sen$.estimate
    spe <- spec(yhatT,class, .pred_class)
    spe <- spe$.estimate
    #rocAUC <- roc_auc(yhatT,class, .pred_class) #not working for some reason. Wants pred_class as numeric?
    
    #add some identifiers
    mod_name <- class(model1)[1]
    sp <- names(X[i])
    
    #save all the metrics
   # mod1_perf[[i]] = c( sp, mod_name, mathews, sen, spe))

  
  list(mod1_k = mod1_k, yhat = yhat, resid = resid) 
    })  

  #cant get performance to work
  #mod1_perf <- do.call(rbind, mod1_perf)
  #colnames( mod1_perf) <- c('Response','Model','MCC', 'Sensitivity', 'Specificity')
  #mod1_perf <- as.data.frame( mod1_perf )
 
}

#---------------------------------------------------------------------------------------------------------------
#yhats <- mrTidyPredict(X=X,Y-Y, model1=model1)
#str(yhats)
#yhats[[1]]$resid
#---------------------------------------------------------------------------------------------------------------
StackPredictions = function(yhats, data.test, covariance_mod = 'gbm'){
  
  # Determine the total number of outcome variables in the yhat data
  n_variables <- length(yhats)
  
  # Check that the covariance_mod specification is allowed
  if(!(covariance_mod %in% c('gbm', 'gam', 'lm')))
    stop('Please select one of the three covariance model options:
         "gbm", "gam", "lm"')
  
  # Need a check for total outcomes so we don't store huge numbers of trees if the 
  # covariance_mod is specified as gbm
  if(n_variables > 50 & covariance_mod == 'gbm'){
    warning('gbm not allowed for data with >50 outputs due to memory constraints. Setting covariance_mod to gam',
            call. = FALSE)
    covariance_mod <- 'gam'
  }
  
  # For each outcome, use model 1 residuals (resid from the yhats list) as outcomes and 
  # predictions from each other outcome's mod 1 (yhat) as predictors to learn the covariance matrix
  cat('Fitting sequential', covariance_mod, 'models to learn the multivariate covariance matrix...\n')
  rhats <- lapply(seq(1, n_variables), function(species){
    outcome <- yhats[[species]]$resid
    predictors <- do.call(cbind, purrr::map(yhats, 'yhat'))[, -species]
    colnames(predictors) <- paste0('pred_', seq(1, ncol(predictors)))
    
    if(covariance_mod == 'gbm'){
    # Use a boosted regression tree to learn complex, nonlinear covariance relationships
      
      #can use the TidyModel syntax here too
    mod2 <- gbm::gbm(Y ~., data = cbind(data.frame(Y = outcome),
                                          as.matrix(predictors)),
                       n.trees = 1000, shrinkage = 0.01, interaction.depth = 3,
                       cv.folds = 10, keep.data = T)
    }
    
    if(covariance_mod == 'gam'){
      # Use a GAM to learn the covariance matrix
      # Include each other outcome as a smooth predictor and use penalisation to
      # force coefficients to zero if support is limited. No interactions here
      max_knots <- apply(predictors, 2, function(x) length(unique(x)))
      
      # Need to set terms with few max_knots to be linear terms
      predictor_formula <- matrix(NA, nrow = 1, ncol = ncol(predictors))
      for(i in 1:ncol(predictors)){
        predictor_formula[1, i] <- ifelse(max_knots[i] < 5, 
                                          colnames(predictors)[i],
                                          paste0("s(", colnames(predictors)[i],", 
                                                bs = 'cs', k = 5",
                                                ")"))
      }
      
      predictor_terms <- paste0(predictor_formula[1,],
                                collapse = '+')
      gam_formula <- as.formula(paste0("Y~", 
                                       predictor_terms))
      mod2 <- mgcv::gam(gam_formula, data = cbind(data.frame(Y = outcome),
                                            as.matrix(predictors)))
    }
    
    if(covariance_mod == 'lm'){
      # Use a simple linear model to learn the covariance matrix
      lm_formula <- as.formula(paste0("Y~", 
                                       paste0(colnames(predictors),collapse = "+")))
      mod2 <- lm(lm_formula, data = cbind(data.frame(Y = outcome),
                                                    as.matrix(predictors)))
    }
    
    # Return the model as well as the 'multivariate residual adjustments'
    list(mod2 = mod2, pred_names = colnames(predictors))
  })
  
  # Now to generate the stacked predictions
  cat('Adjusting original predictions using the proportion of the inverted distance metric...\n')
  combiner_preds <- lapply(seq(1, n_variables), function(species){
    
    # Predict mod 1 for each species using the out-of-sample test data
    mod1_test_preds <- do.call(cbind, lapply(seq(1, n_variables), function(y){
      # may need to find class of mod1 in order to specify the predict function: class(yhats[[1]]$mod1)
      predict(yhats[[y]]$mod1, newdata = data.test, type = 'response') 
    }))
    
    # Predict mod 2 'adjustments' using the test mod 1's predictions
    mod2_test_preds <- do.call(cbind, lapply(seq(1, n_variables), function(y){
      test_preds <- data.frame(mod1_test_preds[,-y])
      colnames(test_preds) <- rhats[[y]]$pred_names
      
      if(covariance_mod == 'gbm'){
        adj_preds <- gbm::predict.gbm(rhats[[y]]$mod2, newdata = test_preds, verbose = F)
      }
      if(covariance_mod == 'gam'){
        adj_preds <- mgcv::predict.gam(rhats[[y]]$mod2, newdata = test_preds)
      }
      if(covariance_mod == 'lm'){
        adj_preds <- predict(rhats[[y]]$mod2, newdata = test_preds)
      }
      adj_preds
    }))
    
    # Adjust predictions using the Xing et al 'proportion of the inverted distance' metric
    # and convert them to the outcome probability scale
    d1 <- abs(1 / (mod2_test_preds[,species] - sqrt((-2 * log(mod1_test_preds[,species])))))
    d0 <- abs(1 / (mod2_test_preds[,species] + sqrt((-2 * log(mod1_test_preds[,species])))))
    final_preds <- d1 / (d1 + d0)
    
    # Return both the probability predictions and the binary predictions
    final_preds_binary <- ifelse(final_preds > 0.4999, 1, 0)
    list(probability_preds = final_preds,
         binary_preds = final_preds_binary)
  })
  
  list(probability_preds = do.call(cbind, purrr::map(combiner_preds, 'probability_preds')),
       binary_preds = do.call(cbind, purrr::map(combiner_preds, 'binary_preds')))
}
#---------------------------------------------------------------------------------------------------------
#test <- StackPredictions(yhats, data.test=data_test, covariance_mod = 'gam')

#get this error:  Error in hardhat::forge(new_data, blueprint) : 
#argument "new_data" is missing, with no default