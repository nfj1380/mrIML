#'Constructs stacked model using the covariance matrix as well as original yhats
#'

#test with rhats$rhat_k should be the correct length.


stacked_preds<- function(rhats, yhats){
  
  n_variables <- length(yhats)
  
  # Now to generate the stacked predictions
  combiner_preds <- lapply(seq(1, n_variables), function(i){
    
    
    mod1_test_preds <- yhats %>%  purrr::map(pluck('yhatT')) #select test data
    
    mod1_test_preds <- do.call(cbind,  mod1_test_preds)
    
    m <- mod1_test_preds[, grep("class", names( mod1_test_preds))] #matching colnames was an issue
    
    # mod1_test_preds<-  m %>% select(!contains("_")) #original data
    mod1_test_preds<-  m %>% select(contains("_"))# %>% as.matrix() %>% as.numeric()#matrix makes it not a factor anymore but may lead to other issues
    mod1_test_preds <- mutate_all( mod1_test_preds  , function(x) as.numeric(as.character(x)))
   
    
    # Predict mod 2 'adjustments' using the test mod 1's predictions
    mod2_test_preds <- do.call(cbind, lapply(seq(1, n_variables), function(y){
      
      test_preds <- as.data.frame(mod1_test_preds[,-y]) ###   
      colnames(test_preds) <- rhats[[y]]$pred_names ###
      test_preds <- mutate_all( test_preds , function(x) as.numeric(as.character(x))) #this is now in the function
    
      #adj_preds <- predict(rhats[[y]]$mod2, new_data = test_preds, verbose = F)
      adj_preds <- predict(rhats[[y]]$mod2, new_data = test_preds) ###
      
      adj_preds
    }))
    
    # Adjust predictions using the Xing et al 'proportion of the inverted distance' metric
    # and convert them to the outcome probability scale
    d1 <- abs(1 / (mod2_test_preds[,i] - sqrt((-2 * log(mod1_test_preds[,i])))))
    d0 <- abs(1 / (mod2_test_preds[,i] + sqrt((-2 * log(mod1_test_preds[,i])))))
    
    
    final_preds <- d1 / (d1 + d0)
    
    # Return both the probability predictions and the binary predictions
    final_preds_binary <- ifelse(final_preds > 0.4999, 1, 0)
    list(probability_preds = final_preds,
         binary_preds = final_preds_binary)
  })
  
 # list(probability_preds = do.call(cbind, purrr::map(combiner_preds, 'probability_preds')),
  #     binary_preds = do.call(cbind, purrr::map(combiner_preds, 'binary_preds')))
}