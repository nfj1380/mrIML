#' Multiple Response Variable Selection
#' 
#' Identify predictors not useful in the prediction of multiple response variables using the Boruta algorithm.
#' 
#' @param X Dataframe of predictors. (Optional)
#' @param X1 Dataframe of additional predictors for predicting the response variables. (Optional)
#' @param Y Dataframe with response variables, where column names are the response variable names.
#' @return A list where each element corresponds to a response variable, containing the Boruta result object for that response variable.
#' @examples
#' X <- data.frame(matrix(rnorm(1000), ncol = 10))
#' X1 <- data.frame(matrix(rnorm(100), ncol = 5))
#' Y <- data.frame(response1 = rnorm(100), response2 = rnorm(100))
#' results <- mrVariableSelect(X, X1, Y)
#' 
mrVariableSelect <- function(X = NULL, X1 = NULL, Y, drop_threshold=0.75) {
  # Check if both X and X1 are provided
  if (!is.null(X) && !is.null(X1)) {
    predictors <- cbind(X, X1)
  } else if (!is.null(X)) {
    predictors <- X
  } else if (!is.null(X1)) {
    predictors <- X1
  } else {
    stop("Either 'X' or 'X1' must be provided.")
  }
  
  # Extract response variable names
  response_vars <- colnames(Y)
  
  response_num <- ncol(Y)
  
  # Initialize an empty list to store results
  results <- list()
  
  # Iterate over response variables
  for (response_var in response_vars) {
    
    # Combine predictors and response variable
    data <- cbind(Y[[response_var]], predictors )
    
    resp_name <- names(Y[response_var])
    
    df1 <- data %>%
      rename(!!resp_name := colnames(data)[1])
    
    # Prepare formula for Boruta
    formula <- as.formula(paste( resp_name , "~ ."))
    
    names_col <- data.frame(rep(resp_name, ncol(X)))
    # Run Boruta algorithm
    boruta_result <- Boruta(formula, data = df1)
    
    # Store results
    results[[response_var]] <- boruta_result
    
  
    # Plot how many times each predictor has been dropped. Reject is more conservative than accept
    results[[response_var]] <- data.frame(predictor = names(boruta_result$finalDecision),
                            dropped = boruta_result$finalDecision == "Rejected", respo =  names_col)
    
  }
  
  reject_data <- do.call( rbind, results) %>% 
    rename(response =  rep.resp_name..ncol.X..)
  
  dropped_prop <-  reject_data %>%
    group_by(predictor) %>%
    summarise(dropped_prop = sum(dropped)/response_num)
  
  # Plot the count of dropped predictors
  p1 <- ggplot(dropped_prop, aes(x = predictor, y = dropped_prop, fill = dropped_prop)) +
    geom_bar(stat = "identity") +
    labs(title = "Predictors dropped by Boruta across responses",
         x = "Predictor",
         y = "Proportion rejected") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  varSelectX <- dropped_prop %>% filter(dropped_prop < drop_threshold) 
  
  # Function to filter columns based on row names
  filter_columns <- function(df, row_names_df) {
    cols_to_keep <- row_names_df$predictor
    df_filtered <- df %>% select(cols_to_keep)
    return(df_filtered)
  }
  
  # Filtering columns of df1 based on row names of df2
  filtered_X <- filter_columns(X, varSelectX)
  

  return(list(p1, filtered_X))
}
