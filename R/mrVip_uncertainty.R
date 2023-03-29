#'Wrapper to estimate and plot model-agnostic variable importance with a measure of uncertainty for multi-response models. 
#'@param Y A \code{dataframe} is response variable data (species, OTUs, SNPs etc).
#'@param X A \code{dataframe} represents predictor or feature data.
#'@param dummy A \code{logical} 'TRUE or FALSE'. 
#'@param Model 1 A \code{list} can be any model from the tidy model package. See examples.
#'@param tune_grid_size A \code{numeric} sets the grid size for hyperparamter tuning. Larger grid sizes increase computational time.
#'@param k A \code{numeric} sets the number of folds in the 10-fold cross-validation. 10 is the default.
#'@seed A \code{numeric} as these models have a stochastic component, a seed is set to make to make the analysis reproducible. Defaults between 100 million and 1.
#'@param mode \code{character}'classification' or 'regression' i.e., is the generative model a regression or classification?

#'@details Calculates variable importance across multiple runs to quantify uncertainty in model estimates.
#'# this is particularly useful for smaller unbalanced data sets where the vip measure of variable
#'#importance can be a bit unstable.
#' @example 
#' # test <- mrVip_mutlirun(X=X,
#'                           Y=Y,
#'                           ity=4, #runs the model 4 times and summarizes. 
#'                          Model=model_rf,
#'                          mode='classification', 
#'                           tune_grid_size=5,
#'                          seed = sample.int(1e8, 1))#'@export 

mrVip_uncertainty<- function (X, Y, Model,
                            mode='regression',ity=5, transformY='log',
                            dummy=FALSE, tune_grid_size= 10, k=10,
                            seed = sample.int(1e8, 1) ) { 

    
    internal_fit_function <- function( i ){
      
    yhats_rf <- mrIMLpredicts(X=X,
                              Y=Y, 
                              Model=model_rf,
                              balance_data='no', 
                              mode='classification', 
                              tune_grid_size=5,
                              seed = sample.int(1e8, 1) ) #this will make sure dif seeds are used
    
    VI <- mrVip(yhats_rf, X=X)   
    
    }
    
    im <- future_lapply(seq(1,ity), internal_fit_function, future.seed = TRUE)

    ImpGlobal <- as.data.frame(do.call(rbind, im))
    
    #make wide to long
    data_vi <- ImpGlobal %>%
      gather(names(X), key = Xvar, value = importance)
    
    #boxplot
    data_vi %>% 
      ggplot(aes(x=reorder(Xvar,importance), y=importance, fill=Xvar)) +
      geom_boxplot()+
      coord_flip()+
      theme_bw()+
      labs(x="Features", y="Importance")+
      scale_fill_brewer(palette="BuPu")#from cbind


}