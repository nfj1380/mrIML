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



mrInteractions <- function(yhats, X, Y, num_bootstrap = 1, feature = feature, top.int = 10) {
  
  internal_fit_function <- function(k) {
    features <- colnames(yhats[[k]]$data)[-1]
    n <- nrow(yhats[[k]]$data)
    int_raw <- list()
    for (i in 1:num_bootstrap) {
      bootstrap_sample <- if (num_bootstrap > 1) {
        yhats[[k]]$data[sample(1:n, replace = TRUE), ]
      } else {
        yhats[[k]]$data
      }
      model_fit <- if (num_bootstrap > 1) {
        wflow <- yhats[[k]]$last_mod_fit %>% tune::extract_workflow()
        wflow$data <- bootstrap_sample
        fit(wflow, data = bootstrap_sample)
      } else {
        yhats[[k]]$mod1_k %>% tune::extract_fit_parsnip()
      }
      
      pred_fun <- function(m, dat) {
        predict(m, dat[, colnames(bootstrap_sample)[-1], drop = FALSE], type = "prob")$.pred_1
      }
      s <- hstats(model_fit, v = names(yhats[[k]]$data_train)[-1], 
                  X = yhats[[k]]$data_train, pred_fun = pred_fun, 
                  n_max = 300, pairwise_m = length(names(yhats[[k]]$data_train)[-1]), 
                  threeway_m = 0, verbose = FALSE)
      overall <- data.frame(response = names(Y[k]), overall = h2(s)[[1]], bs = i)
      one_way <- data.frame(one_way = h2_overall(s, plot = FALSE)[[1]]) %>% 
        rownames_to_column("predictor")
      metadata <- data.frame(response = rep(names(Y[k]), nrow(one_way)), bstrap = rep(i, nrow(one_way)))
      one_way_df <- cbind(one_way, metadata)
      two_way <- data.frame(two_way_int = h2_pairwise(s, plot = FALSE)[[1]]) %>% 
        rownames_to_column("predictor")
      meta_data2 <- data.frame(response = rep(names(Y[k])), bstrap = rep(i, nrow(two_way)))
      two_way_df <- cbind(two_way, meta_data2)
      interaction_objects <- list(overall_int = overall, 
                                  one_way_int = one_way_df, two_way_int = two_way_df)
      int_raw[[i]] <- interaction_objects
    }
    return(int_raw)
  }
  
  bstraps_int_list <- future_lapply(seq_along(yhats), internal_fit_function, future.seed = TRUE)
  
  # Combine bootstrap results
  overall_int_final <- do.call(rbind, lapply(bstraps_int_list, function(sublist) {
    extracted_overall_int <- map(sublist, pluck, "overall_int")
    do.call(rbind, extracted_overall_int)
  }))
  
  # Plot overall interactions
  p1 <- ggplot(overall_int_final, aes(x = reorder(response, -overall), y = overall)) + 
    geom_boxplot() + labs(title = "Overall interactions", x = "Response", y = "Overall") + theme_bw()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Combine one-way interaction results
  overall_one_way_final <- do.call(rbind, lapply(bstraps_int_list, function(sublist) {
    extracted_one_way_int <- map(sublist, pluck, "one_way_int")
    do.call(rbind, extracted_one_way_int)
  }))
  
  # Filter and plot one-way interactions
  filtered_one_way <- overall_one_way_final %>% filter(response == feature) %>% 
    slice_max(order_by = one_way, n = top.int)
  p2 <- ggplot(filtered_one_way, aes(x = reorder(predictor, -one_way), y = one_way)) + 
    geom_boxplot() + labs(title = paste(feature, "one-way interactions", sep = " "), 
                          x = feature, y = paste(feature, "interaction importance", sep = " ")) + theme_bw()
  
  # Combine two-way interaction results
  overall_two_way_final <- do.call(rbind, lapply(bstraps_int_list, function(sublist) {
    extracted_two_way_int <- map(sublist, pluck, "two_way_int")
    do.call(rbind, extracted_two_way_int)
  }))
  
  # Filter and plot two-way interactions
  filtered_two_way <- overall_two_way_final %>% filter(response == feature) %>% 
    slice_max(order_by = two_way_int, n = top.int)
  
  p3 <- ggplot(filtered_two_way, aes(x = reorder(predictor, -two_way_int), y = two_way_int)) + 
    geom_boxplot() + labs(title = paste(feature, "two-way interactions", sep = " "), 
                          x = feature, y = paste(feature, "interaction importance", sep = " ")) + 
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(list(p1, p2, p3, overall_int_final, overall_one_way_final, overall_two_way_final))
}
