#' Calculate and visualize feature interactions
#'
#' This function calculates and visualizes interactions in the model using bootstrapping.
#' It provides overall, one-way, and two-way interactions for specified features.
#'
#' @param yhats A list of model predictions.
#' @param X The predictor data.
#' @param Y The response data.
#' @param num_bootstrap The number of bootstrap samples to generate (default: 1).
#' @param feature The feature for which interactions need to be calculated.
#' @param top.int The number of top interactions to display (default: 10).
#'
#' @return A list containing the visualizations for overall, one-way, and two-way interactions, as well as the interaction dataframes.
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' #set up analysis
#' Y <- dplyr::select(Bird.parasites, -scale.prop.zos)%>% 
#' dplyr::select(sort(names(.)))#response variables eg. SNPs, pathogens, species....
#' X <- dplyr::select(Bird.parasites, scale.prop.zos) # feature set

#' X1 <- Y %>%
#' dplyr::select(sort(names(.)))
#'model_rf <- 
#' rand_forest(trees = 100, mode = "classification", mtry = tune(), min_n = tune()) %>% #100 trees are set for brevity. Aim to start with 1000
#' set_engine("randomForest")
#' yhats_rf <- mrIMLpredicts(X=X, Y=Y,
#'X1=X1,'Model=model_rf , 
#'balance_data='no',mode='classification',
#'tune_grid_size=5,seed = sample.int(1e8, 1),'morans=F,
#'prop=0.7, k=5, racing=T) 
#'int_ <- mrInteractions(yhats=yhats_rf, X, Y, num_bootstrap=10,
#'feature = 'Microfilaria', top.int=10)
#'int_[[1]] # overall plot
#'int_[[2]] # individual plot for the response of choice 
#'int_[[3]] #two way plot# }


mrInteractions <- function(yhats, X, Y, num_bootstrap = 1,  
                           feature = feature, top.int=10) {
  
  n_response <- length(yhats)
  
  pb <- txtProgressBar(min = 0, max = n_response, style = 3) 
  
  internal_fit_function <- function(k) {
    setTxtProgressBar(pb, k)
    
    features <- colnames(yhats[[k]]$data)[-1]
    n <- nrow(yhats[[k]]$data)
    
    int_raw <- list()
    
    if (is.na(num_bootstrap) || num_bootstrap == 1 ) {
      bootstrap_sample <- yhats[[k]]$data
    }
    
    for (i in 1:num_bootstrap) {
      if (num_bootstrap > 1) {
        bootstrap_sample <- yhats[[k]]$data[sample(1:n, replace = TRUE), ]
        wflow <- yhats[[k]]$last_mod_fit %>% tune::extract_workflow()
        wflow$data <- bootstrap_sample
        model_fit <- fit(wflow, data = bootstrap_sample)
      } else {
        model_fit <- yhats[[k]]$mod1_k %>% extract_fit_parsnip()
      }
      
      metrics <- list(
        logloss = MetricsWeighted::logLoss,
        `ROC AUC` = MetricsWeighted::AUC,
        `% Dev Red` = MetricsWeighted::r_squared_bernoulli
      )
      
      var_names <- names(yhats[[k]]$data)[-1]
      
      pred_fun <- function(m, dat) {
        predict(
          m, dat[, colnames(bootstrap_sample)[-1], drop = FALSE],
          type = "prob"
        )$`.pred_1`
      }
      
      s <- hstats(model_fit, v = names(yhats[[k]]$data_train)[-1],
                  X = yhats[[k]]$data_train, pred_fun = pred_fun, n_max = 300, 
                  pairwise_m = length(names(yhats[[k]]$data_train)[-1]), threeway_m = 0, verbose=F )
      
      overall <- data.frame(response = names(Y[k]), overall = h2(s), bs = i)
      
      one_way <- as.data.frame(h2_overall(s, plot = FALSE)) %>% 
        rownames_to_column('predictor')
      
      metadata <-  data.frame(response = rep(names(Y[k]), nrow(one_way)), bstrap = rep(i, nrow(one_way)))
      
      one_way_df <- cbind(one_way, metadata)
      
      two_way <- as.data.frame(h2_pairwise(s, plot = FALSE), top_m = inf) %>% 
        rownames_to_column('predictor')
      
      meta_data2 <-  data.frame(response = rep(names(Y[k])), bstrap = rep(i, nrow(two_way)))
      
      two_way_df <- cbind(two_way, meta_data2)
      
      interaction_objects <- list(
        overall_int = overall,
        one_way_int = one_way_df,
        two_way_int = two_way_df
      )
      
      int_raw[[i]] <- interaction_objects
    }
    
    return(int_raw)
  }
  
  bstraps_int_list <- future_lapply(seq(1, n_response), internal_fit_function, future.seed = TRUE)
  
  ######################  
  # overall interactions
  ######################
  
  #how much variability do interactions account for?
  bs_list_overall <- lapply(bstraps_int_list, function(sublist) {
    extracted_overall_int <- map(sublist, pluck, "overall_int")
    do.call(rbind, extracted_overall_int)
  })
  
  overall_int_final <- do.call(rbind, bs_list_overall)
  
  #showing just the top interactions only
  top_int_overall <-  overall_int_final %>%
    group_by(response) %>%
    summarize(avg_Int = mean(overall), .groups = "drop")
  
  top_int_overall_ordered <- top_int_overall %>% 
    arrange(desc(avg_Int)) %>%  # Sort by in descending order
    dplyr::slice(1:top.int) 
  top_names <-  top_int_overall_ordered$response
  
  overall_int_final_top <- overall_int_final %>% 
    dplyr::filter(response %in% top_names)
  
  p1 <- ggplot(overall_int_final_top, aes(x = reorder(response, -overall), y = overall)) +
    geom_boxplot() +
    labs(title = "Overall interactions", x = "Response", y = "Overall") +
    theme_bw()
  
  ######################  
  #one-way interactions
  ######################
  
  bs_list_one_way <- lapply(bstraps_int_list, function(sublist) {
    extracted_one_way_int <- map(sublist, pluck, "one_way_int")
    do.call(rbind, extracted_one_way_int)
  })
  
  overall_one_way_final <- do.call(rbind, bs_list_one_way)
  
  filtered_one_way <- overall_one_way_final %>% 
    dplyr::filter(response == feature) 
  
  top_int_one_way <-   filtered_one_way  %>%
    group_by(response,predictor) %>%
    summarize(avg_Int = mean(V1), .groups = "drop")
  
  top_int_one_way_ordered <- top_int_one_way %>% 
    arrange(desc(avg_Int)) %>%  # Sort by in descending order
    dplyr::slice(1:top.int) 
  
 top_names_one_way <-  top_int_one_way_ordered$predictor
  
 one_way_int_final_top <- filtered_one_way %>% ###made change here
   dplyr::filter(predictor %in% top_names_one_way)
 
 p2 <- ggplot(one_way_int_final_top, aes(x = reorder(predictor, -V1), y = V1)) +
   geom_boxplot() +
   labs(title = paste(feature,"one-way interactions", sep= " "), x = feature,
        y = paste(feature, "interaction importance", sep = " ")) +
   theme_bw()
 
 
  #community. Takes the average across all bstraps
  avg_v_by_response <- overall_one_way_final %>%
    group_by(response, predictor) %>%
    summarize(avg_V1 = mean(V1), .groups = "drop") %>% 
    arrange(desc(avg_V1)) %>% 
    dplyr::slice(1:top.int) 
  
  p2_com <- ggplot(avg_v_by_response, aes(x = reorder(predictor, -avg_V1), y = avg_V1)) +
    geom_boxplot() +
    labs(title = "Community-level one-way interactions", x = 'Predictors',
         y ="Community interaction importance") +
    theme_bw()
  
  combined_plot_one_way <- plot_grid(p2, p2_com, ncol = 2)

######################  
#two-way interactions
######################
  
bs_list_two_way <- lapply(bstraps_int_list, function(sublist) {
  extracted_two_way_int <- map(sublist, pluck, "two_way_int")
  do.call(rbind, extracted_two_way_int)
})

overall_two_way_final <- do.call(rbind, bs_list_two_way)
  
filtered_two_way <- overall_two_way_final %>% 
  dplyr::filter(response == feature) 

top_int_two_way <-   filtered_two_way  %>%
  group_by(response,predictor) %>%
  summarize(avg_Int = mean(V1), .groups = "drop")

top_int_two_way_ordered <- top_int_two_way %>% 
  arrange(desc(avg_Int)) %>%  # Sort by in descending order
  dplyr::slice(1:top.int) 
top_names_two_way <-  top_int_two_way_ordered$predictor

two_way_int_final_top <- filtered_two_way %>% 
  filter(predictor %in% top_names_two_way )

p3 <- ggplot(two_way_int_final_top, aes(x = reorder(predictor, -V1), y = V1)) +
  geom_boxplot() +
  labs(title = paste(feature,"two-way interactions", sep= " "), x = feature,
       y = paste(feature, "interaction importance", sep = " ")) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))

avg_v_by_response_two_way <- overall_two_way_final %>%
  group_by(response, predictor) %>%
  summarize(avg_V1 = mean(V1), .groups = "drop") %>%  # Sort by in descending order
  dplyr::slice(1:top.int) 

p3_com <- ggplot(avg_v_by_response_two_way, aes(x = reorder(predictor, -avg_V1), y = avg_V1)) +
  geom_boxplot() +
  labs(title = "Community-level two-way interactions", x = 'Predictors',
       y = "Community interaction importance") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
  
combined_plot_two_way <- plot_grid(p3, p3_com, ncol = 2)

  return(list(p1, combined_plot_one_way, combined_plot_two_way,
              overall_int_final, overall_one_way_final,
              overall_two_way_final))
}
