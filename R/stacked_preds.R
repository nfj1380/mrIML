#'Constructs stacked model using the covariance matrix as well as original yhats
#'

#test with rhats$rhat_k should be the correct length.


stacked_preds<- function(rhats, yhats){
  
  # Now to generate the stacked predictions
  combiner_preds <- lapply(seq(1, 4), function(species){
    
    
    mod1_test_preds <- yhats %>%  purrr::map(pluck('yhatT')) #select test data
    
    mod1_test_preds <- do.call(cbind,  mod1_test_preds)
    
    m <- mod1_test_preds[, grep("class", names( mod1_test_preds))] #matching colnames was an issue
    
    # mod1_test_preds<-  m %>% select(!contains("_")) #original data
    mod1_test_preds<-  m %>% select(contains("_")) %>% as.matrix() #matrix makes it not a factor anymore but may lead to other issues
    
    # Predict mod 2 'adjustments' using the test mod 1's predictions
    mod2_test_preds <- do.call(cbind, lapply(seq(1, n_variables), function(y){
      test_preds <- mod1_test_preds[,-y] #not a dataframe anymore  - could be an issue. Otherwise is was a factor
      colnames(test_preds) <- rhats[[y]]$pred_names
      
      #adj_preds <- predict(rhats[[y]]$mod2, new_data = test_preds, verbose = F)
      adj_preds <- predict(rhats[[y]]$mod2, new_data = test_preds) 
      
      adj_preds
    }))
    
    # Adjust predictions using the Xing et al 'proportion of the inverted distance' metric
    # and convert them to the outcome probability scale
    d1 <- abs(1 / (mod2_test_preds[,species] - sqrt((-2 * log(mod1_test_preds[,species]))))) #something isn't quite right here. mod1_test_preds[,species]) is a vector
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