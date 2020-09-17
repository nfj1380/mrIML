#Function to stack the response models and compare accuracy of stacked versus unstacked models

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

  #--------------------------------------------------------------------------------------------------------- 
 
  
  # Predict mod 1 for each species using the out-of-sample test data
   mod1_test_preds <- do.call(cbind, lapply(seq(1, 4), function(y){
    predict(yhats[[1]]$mod1_k, newdata = test.dat, type = 'response') 
  }))
  
  # Predict mod 2 'adjustments' using the test mod 1's predictions
  mod2_test_preds <- do.call(cbind, lapply(seq(1, 4), function(y){
    test_preds <- data.frame(mod1_test_preds[,-y])
    colnames(test_preds) <- rhats[[y]]$pred_names
    predict(rhats[[y]]$mod2_k, newdata = test_preds)
  }))
  
  # Adjust predictions using the 'proportion of the inverted distance' metric
  # and convert them to the probability scale
  d1 <- abs(1 / (mod2_test_preds[,species] - sqrt((-2 * log(mod1_test_preds[,species])))))
  d0 <- abs(1 / (mod2_test_preds[,species] + sqrt((-2 * log(mod1_test_preds[,species])))))
  
  final_preds <- d1 / (d1 + d0)
  final_preds <- ifelse(final_preds > 0.4999, 1, 0)
  final_preds
})

# Calculate prediction accuracies for the stacked model and compare to the 
# simpler univariate model (i.e. model 1 only) for each species
lapply(seq(1, 4), function(species){
  stacked_stats <- caret::confusionMatrix(factor(combiner_preds[[species]]), 
                                          factor(as.matrix(test.dat)[, species]))
  
  single_preds <- predict(yhats[[species]]$mod1_k, newdata = test.dat, type = 'response')
  single_preds <- ifelse(single_preds > 0.4999, 1, 0)
  single_stats <- caret::confusionMatrix(factor(single_preds), 
                                         factor(as.matrix(test.dat)[, species]))
  list(stacked_stats = round(stacked_stats$overall, 4),
       single_stats = round(single_stats$overall, 4))
})
}
# Out-of-sample prediction accuracy improves for nearly all species