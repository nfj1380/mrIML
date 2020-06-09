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
#'data("Bird.parasites")
#'
#'# A simple split of the data into training and testing subsets
#'train.dat <- Bird.parasites[1:350, ]
#'test.dat <- Bird.parasites[351:449, ]
#'
#'# Run model 1 for each parasite; a simple logistic regression with a single covariate
#'# in this case but note that model 1 can be any model of the user's choice, 
#'# from simple regressions to complex hierarchical or deep learning models.
#'
#'# Different structures can also be used for each species to handle mixed multivariate outcomes
#'yhats <- lapply(colnames(Bird.parasites)[1:4], function(species){
#'
#'# Fit model one for each parasite; can easily modify this so that the user
#'# can specify the formula necessary for each species as a list of formulas
#'
#'mod1 <- glm(formula(paste0(species,' ~ scale.prop.zos')), 
#'family = 'binomial', data = train.dat)
#'
#'# Calculate probability predictions for the fitted training data
#'yhat <- predict(mod1, type = 'response')
#'
#'# Calculate deviance residuals
#'# To calc residuals consistently for gbm or others for binomial models we need the formula
#'# train.dat[, species] equals binary observations; yhat equals model probability predictions
#'resid <- ifelse(train.dat[, species] == 1, 
#'sqrt((-2 * log(yhat))), 
#'-1 * (sqrt((-2*log(1 - yhat)))))
#'resid <- as.vector(resid)
#'#resid <- residuals(mod1, "deviance")
#'list(mod1 = mod1, yhat = yhat, resid = resid)
#'})
#'
#'# A GAM example for the covariance matrix
#'gam_stack <- StackPredictions(yhats, data.test = test.dat, covariance_mod = 'gam')
#'
#'# Calculate prediction accuracies for the stacked model and compare to the 
#'# simpler univariate model (i.e. model 1 only) for each species
#'lapply(seq(1, 4), function(species){
#'stacked_stats <- caret::confusionMatrix(factor(gam_stack$binary_preds[, species]), 
#'factor(as.matrix(test.dat)[, species]))
#'single_preds <- predict(yhats[[species]]$mod1, newdata = test.dat, type = 'response')
#'single_preds <- ifelse(single_preds > 0.4999, 1, 0)
#'single_stats <- caret::confusionMatrix(factor(single_preds), 
#'factor(as.matrix(test.dat)[, species]))
#'list(stacked_stats = round(stacked_stats$overall, 4),
#'single_stats = round(single_stats$overall, 4))
#'})
#'
#'@export
#'
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
    mod2_k <- gbm::gbm(Y ~., data = cbind(data.frame(Y = outcome),
                                          as.matrix(predictors)),
                       n.trees = 1000, shrinkage = 0.01, interaction.depth = 3,
                       cv.folds = 10, keep.data = T)
    rhat_k <- predict(mod2_k)
    }
    
    if(covariance_mod == 'gam'){
      # Use a GAM to learn the covariance matrix
      # Include each other outcome as a smooth predictor and use penalisation to
      # force coefficients to zero if support is limited. No interactions here
      gam_formula <- as.formula(paste0("Y~", 
                                       paste0("s(",colnames(predictors),
                                              ", bs = 'cs', k = 5)",collapse = "+")))
      mod2_k <- mgcv::gam(gam_formula, data = cbind(data.frame(Y = outcome),
                                            as.matrix(predictors)))
      rhat_k <- mgcv::predict.gam(mod2_k, outcome = 'response')
    }
    
    if(covariance_mod == 'lm'){
      # Use a simple linear model to learn the covariance matrix
      lm_formula <- as.formula(paste0("Y~", 
                                       paste0(colnames(predictors),collapse = "+")))
      mod2_k <- lm(lm_formula, data = cbind(data.frame(Y = outcome),
                                                    as.matrix(predictors)))
      rhat_k <- predict(mod2_k, outcome = 'response')
    }
    
    # Return the model as well as the 'multivariate residual adjustments'
    list(rhat_k = rhat_k,
         mod2_k = mod2_k, pred_names = colnames(predictors))
  })
  
  # Now to generate the stacked predictions
  cat('Adjusting original predictions using the proportion of the inverted distance metric...\n')
  combiner_preds <- lapply(seq(1, n_variables), function(species){
    
    # Predict mod 1 for each species using the out-of-sample test data
    mod1_test_preds <- do.call(cbind, lapply(seq(1, n_variables), function(y){
      # may need to find the class in order to specify the predict function class(yhats[[1]]$mod1)
      predict(yhats[[y]]$mod1, newdata = data.test, type = 'response') 
    }))
    
    # Predict mod 2 'adjustments' using the test mod 1's predictions
    mod2_test_preds <- do.call(cbind, lapply(seq(1, 4), function(y){
      test_preds <- data.frame(mod1_test_preds[,-y])
      colnames(test_preds) <- rhats[[y]]$pred_names
      
      if(covariance_mod == 'gbm'){
        adj_preds <- gbm::predict.gbm(rhats[[y]]$mod2_k, newdata = test_preds, verbose = F)
      }
      if(covariance_mod == 'gam'){
        adj_preds <- mgcv::predict.gam(rhats[[y]]$mod2_k, newdata = test_preds)
      }
      if(covariance_mod == 'lm'){
        adj_preds <- predict(rhats[[y]]$mod2_k, newdata = test_preds)
      }
      adj_preds
    }))
    
    # Adjust predictions using the 'proportion of the inverted distance' metric
    # and convert them to the outcome scale
    d1 <- abs(1 / (mod2_test_preds[,species] - sqrt((-2 * log(mod1_test_preds[,species])))))
    d0 <- abs(1 / (mod2_test_preds[,species] + sqrt((-2 * log(mod1_test_preds[,species])))))
    
    final_preds <- d1 / (d1 + d0)
    final_preds_binary <- ifelse(final_preds > 0.4999, 1, 0)
    list(probability_preds = final_preds,
         binary_preds = final_preds_binary)
  })
  
  list(probability_preds = do.call(cbind, purrr::map(combiner_preds, 'probability_preds')),
       binary_preds = do.call(cbind, purrr::map(combiner_preds, 'binary_preds')))
}

