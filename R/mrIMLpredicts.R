#'Wrapper to generate multi-response predictive models.
#'@param Y A \code{dataframe} is response variable data (species, OTUs, SNPs etc).
#'@param X A \code{dataframe} represents predictor or feature data.
#'#'@param X1 A \code{dataframe} extra predictor set used in each model. For the MrIML Joint species distribution model (JSDM) this is just a copy of the response data.
#'@param balance_data A \code{character} 'up', 'down' or 'no'. 
#'@param dummy A \code{logical} 'TRUE or FALSE'. 
#'@param Model 1 A \code{list} can be any model from the tidy model package. See examples.
#'@param tune_grid_size A \code{numeric} sets the grid size for hyperparameter tuning. Larger grid sizes increase computational time. Ignored if racing=T.
#'@param racing \code{logical} 'TRUE or FALSE'. If 'TRUE' MrIML performs the grid search using the 'racing' ANOVA method. See https://finetune.tidymodels.org/reference/tune_race_anova.html
#'@param k A \code{numeric} sets the number of folds in the 10-fold cross-validation. 10 is the default.
#'@seed A \code{numeric} as these models have a stochastic component, a seed is set to make to make the analysis reproducible. Defaults between 100 million and 1.
#'@param mode \code{character}'classification' or 'regression' i.e., is the generative model a regression or classification?
#'@param morans \code{logical} 'TRUE or FALSE'. If 'TRUE' global Morans I is calculated for each response

#'@details This function produces yhats that used in all subsequent functions.
#' This function fits separate classification/regression models for each response variable in a data set.  Rows in X (features) have the same id (host/site/population)
#'  as Y. Class imbalance can be a real issue for classification analyses. Class imbalance can be addressed for each
#' response variable using 'up' (upsampling using ROSE bootstrapping), 'down' (downsampling) 
#'or 'no' (no balancing of classes).
#' @example 
#' all_cores <- parallel::detectCores(logical = FALSE)

#'cl <- makePSOCKcluster(all_cores)
#'registerDoParallel(cl)
#'
#' model1 <- 
#' rand_forest(trees = 100, mode = "classification") %>% #this should cope with multinomial data alreadf
#'   set_engine("ranger", importance = c("impurity","impurity_corrected")) %>% #model is not tuned to increase computational speed
#'  set_mode("classification")
#'  
#' yhats <- mrIMLpredicts(X= enviro_variables,Y=response_data, model1=model1, balance_data='no', model='classification',  
#'tune_grid_size=5, k=10, seed = sample.int(1e8, 1)))
#'@export


mrIMLpredicts<- function(X, X1=NULL, Y, Model, balance_data ='no', mode='regression',dummy=FALSE,
                         prop=0.5, morans=F, tune_grid_size= 10, k=10, racing=T, seed = sample.int(1e8, 1) ) { 
  
  n_response<- length(Y)
 
  mod1_perf <- NULL #place to save performance matrix
  
  pb <- txtProgressBar(min = 0, max = n_response, style = 3)
  
  internal_fit_function <- function( i ){
    setTxtProgressBar(pb, i)
    
    if (!is.null(X1)) {
      if (!is.null(X)) {
        data <- as.data.frame(cbind(Y[i], X, X1[-i]), stringsAsFactors = FALSE) ###
      } else {
        data <- as.data.frame(cbind(Y[i], X1[-i]), stringsAsFactors = FALSE)
      }
    } else {
      if (!is.null(X)) {
        data <- as.data.frame(cbind(Y[i], X))
      } else {
        stop("At least one of X or X1 must be specified.")
      }
    }
    
    
    colnames(data)[1] <- c('class') #define response variable for either regression or classification
    
    
    if (mode=='classification'){
      
    data$class<- as.factor(data$class)
    }
  
  
    set.seed(seed)
    
    data_split <- initial_split(data, prop = prop) 
    
    # extract training and testing sets
    data_train <- training(data_split)
    data_test <- testing(data_split)

    #n fold cross validation
    data_cv <- vfold_cv(data_train, v= k) 
    
    if(balance_data == 'down'){ 
      
      data_recipe <- training(data_split) %>%
        recipe(class ~., data= data_train) %>% 
        themis::step_downsample(class) 
    }
    
    if(balance_data == 'up'){
      
      data_recipe <- training(data_split) %>%
        recipe(class ~., data= data_train) %>%
        themis::step_rose(class) #ROSE works better on smaller data sets. SMOTE is an option too.
    }
    
    if(balance_data == 'no'){ 
      
      data_recipe <- training(data_split) %>% 
        recipe(class ~., data= data_train)
    }
    
    if ( dummy == TRUE){
      
      data_recipe <- data_recipe %>% step_dummy(all_nominal(), -all_outcomes(), one_hot = T) #adds dummy variables if needed to any feature that is a factor
    }
    # else{
    #   data_recipe_final <- data_recipe
   # }
    #optional recipe ingredients can be easily added to.
    #step_corr(all_predictors()) %>% # removes all correlated features
    #step_center(all_predictors(), -all_outcomes()) %>% #center features
    #step_scale(all_predictors(), -all_outcomes()) %>% #scale features
    
    mod_workflow <- workflow() %>%
      # add the recipe
      add_recipe(data_recipe) %>%
      # add the model
      add_model(Model)
    
  
    if (racing == TRUE) {
      
      tune_m <- finetune::tune_race_anova(mod_workflow, resamples = data_cv)
      
    } else {
      
      tune_m <- tune::tune_grid(mod_workflow, resamples = data_cv, grid = tune_grid_size)
    }
  
  if (mode=='classification'){
    
    # select the best model
    best_m <- tune_m %>%
      select_best("roc_auc") #try mcc or roc_auc
  }
  
  if (mode=='regression'){
    
    # select the best model
    best_m <- tune_m %>%
      select_best("rmse")
    
  }
  
    
  # final model specification
  final_model <- finalize_workflow(mod_workflow,
                                   best_m )
  # now to fit the model
    mod1_k <- final_model %>%
      fit(data = data_train)
    
    
    # make predictions and calculate deviance residuals. 
    
    if (mode=='classification'){
      
      #predictions
      yhatO <- predict(mod1_k, new_data = data, type='prob' )
      
      yhat <- yhatO$.pred_1
      
      #predictions based on testing data
      yhatT <- predict(mod1_k, new_data = data_test, type='class' ) %>%
        bind_cols(data_test %>% dplyr::select(class))
      
      truth <- as.numeric(as.character(data$class))
      
      #pred_truth <- cbind(yhat, truth)
    
    deviance <- sapply(seq_along(truth), function(j) {
      resid <- ifelse(truth[j] == 1, 
                      sqrt((-2 * log(yhat[j]))), 
                      -1 * (sqrt((-2*log(1 - yhat[j])))))
      return(resid)
    })
    
    deviance_morans <- deviance
    deviance_morans[is.infinite(deviance_morans)] <- 2 #cant have Inf values. It means that the model isnt fitting this
    #data point very well. This is a temporary fix
    }
    
    if (mode=='regression'){
      
      yhatO <- predict(mod1_k, new_data = data_train ) 
      
      yhat <- yhatO$.pred
      
      #predictions based on testing data
      yhatT <- predict(mod1_k, new_data = data_test) # %>% 
      #bind_cols(data_test %>% select(class))
      
      deviance=NULL
      
    }
    
    moran_p <- NA  # Initialize moran_p with a default value
    moran_stat <- NA
    
      if(morans==TRUE){
      
    combined_data <- cbind(spatial_data, deviance_morans)
    
    # Convert data to a SpatialPointsDataFrame or SpatialPolygonsDataFrame
    sp::coordinates(combined_data) <- c("longitude", "latitude")
    
    # Calculate the spatial weights matrix
    spatial_weights <- spdep::dnearneigh(combined_data, d1 = 0, d2 = 1)
    
    #create network
    spatial_weights_listw <- spdep::nb2listw(spatial_weights, style = "W")
    
    # Calculate Moran's I on residuals
    moran_residuals <- spdep::moran.mc(dev, listw =  spatial_weights_listw, nsim = 999)
    
    moran_p <- data.frame(response=names(Y[i]), morans_p=moran_residuals$p.value)
    
    moran_stat <- data.frame(response=names(Y[i]), morans_stat=moran_residuals$stat)
    
    }
    
    
      
    # the last fit. Useful for some functionality

    last_mod_fit <- 
      final_model %>% 
      last_fit(data_split)
    
    #save data
    list(mod1_k = mod1_k, last_mod_fit=last_mod_fit,tune_m=tune_m, data=data, 
         data_testa=data_test, data_train=data_train, yhat = yhat, yhatT = yhatT,
         deviance = deviance, moran_p=moran_p,  moran_stat = moran_stat)
    
    
}

     
    yhats <- future_lapply(seq(1,n_response), internal_fit_function, future.seed = TRUE)
    
   
}


  


