#'Calculates and helps interpret variable importance for mrIML models.
#' @param yhats A list of model predictions.
#' @param mrBootstrap_obj The object containing model bootstrapping results.
#' @param ModelPerf A list containing model performance metrics for each response variable.
#' @param X The predictor data.
#' @param Y The response data.
#' @param X1 A \code{dataframe} extra predictor set used in each model. For the MrIML Joint species distribution model (JSDM) this is just a copy of the response data.
#' @param threshold The threshold for model performance (AUC) below which variables are filtered (default: 0.1).
#' @param global_top_var The number of top global variables to display (default: 10).
#' @param local_top_var The number of top local variables for each response to display (default: 5).
#' @param mode \code{character}'classification' or 'regression' i.e., is the generative model a regression or classification?
#' @return A list containing variable importance data and a combined plot.
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
#'prop=0.7, k=5, racing=T) #
#'ModelPerf <- mrIMLperformance(yhats_rf, Model=model_rf, Y=Y, mode='classification')
#'
#'bs_analysis <- mrBootstrap(yhats=yhats_rf,Y=Y, num_bootstrap = 5)
#'bs_impVI <- mrvip(yhats=yhats_rf, mrBootstrap_obj=bs_analysis, ModelPerf=ModelPerf, 
#'threshold=0.8,  X=X, Y=Y, X1=X1, global_top_var=10,
#'local_top_var=5)
#'bs_impVI[[1]] #importance data bootstrap
#'bs_impVI[[2]] #importance data point estimate
#'bs_impVI[[3]] #importance plot
#'bs_impVI[[4]] #PCA
#'bs_impVI[[5]] #outliers
#'bs_impVI[[6]] #pc loadings
#'bs_impVI[[7]] #PC scores
#'}
#' @export

mrvip <- function(yhats = NULL, mrBootstrap_obj = NULL,  X=X, X1=NULL, Y=Y,
                  mode = 'classification', threshold = 0.1, global_top_var = 10,
                  local_top_var = 5, taxa=NULL, ModelPerf=ModelPerf) {
  

    # If bootstrap_obj is NULL, treat it as a predictive model (yhats)
    
    if (is.null(mode)) {
      stop("Mode argument must be provided when using yhats option.")
    }
    
    #internal_fit_function <- function(i) {
    
    options(dplyr.summarise.inform = FALSE)
    
    if (mode == 'classification') {
      
      pred_fun <- function(m, dat) {
        
        predict(
          m, dat[, colnames(yhats[[i]]$data[-1]), drop = FALSE],
          type = "prob"
        )$`.pred_1`
      }
      
      metrics <- list(
        logloss = MetricsWeighted::logLoss,
        `ROC AUC` = MetricsWeighted::AUC,
        `% Dev Red` = MetricsWeighted::r_squared_bernoulli
      )
      
    } else {
      
      pred_fun <- function(m, dat) {
        predict(m, dat[, colnames(X), drop = FALSE])[[".pred"]]
        
      }
      
      metrics=list(rmse = MetricsWeighted::rmse)
    }
    
    models <- lapply(yhats, `[[`, "mod1_k")
    
    fl_list <- vector("list", length(yhats))
    
    vi_df <- list()  # Initialize a list to store the final data frames
    
    for (i in seq_along(fl_list)) {
      
      fl_list[[i]] <- flashlight(
        
        model = yhats[[i]]$mod1_k,
        label = colnames(Y)[i],
        y = 'class', 
        data = yhats[[i]]$data,
        predict_function = pred_fun,
        metrics = metrics
      )
      
      # Define variable names based on whether X1 is NULL
      if (is.null(X1)) {
        var_names <- colnames(X)
      } else {
        var_names <- colnames(cbind(X, Y[-i]))  
      }
      
      vi_list <- list()  # Initialize a list to store variable importance data frames
      
      # Calculate variable importance for each variable
      for (j in seq_along(var_names)) {
        
        pd_ <- flashlight::light_profile(fl_list[[i]], v = var_names[j])
        
        # Calculate standard deviation of each variable
        sd_value <- sd(pd_$data$value)
        
        # Create a data frame with variable name and standard deviation
        vi_list[[j]] <- data.frame(var = var_names[j], sd_value = sd_value, response = colnames(Y)[i])
      }
      
      # Combine variable importance data frames into one data frame
      vi_df[[i]] <- do.call(rbind, vi_list)
    }
    
    # Combine all variable importance data frames into one table
    vi_table <- do.call(rbind, vi_df)
    
      # Get top variables
      G_target_data_avg <- vi_table %>% 
        dplyr::group_by(var) %>% 
        dplyr::summarise(mean_imp = mean(sd_value)) 
      
      #filter variables
      G_top_vars <- head(G_target_data_avg[order(-G_target_data_avg$mean_imp), ], global_top_var)
      G_target_data_final_df <- vi_table %>% 
        dplyr::filter(var %in% G_top_vars$var)
      
      #first plot
      p1 <- ggplot(G_target_data_final_df, aes(y = reorder(var, sd_value), x = sd_value))+
        geom_col() +
        theme_bw()+
        labs(x = "Importance", y = "Features")+
        theme(
          axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8)
        )
      
      #plot individual taxa 
      
      plot_list <- list()
      
      unique_responses <- unique(G_target_data_final_df$response)
      
      #reduces number of plots - focuses on responses best predicted
      
      if (mode=='classification'){
        
      perf_filter <- ModelPerf[[1]] %>% filter (roc_AUC > threshold)
      
      } else {
        
      perf_filter <- ModelPerf[[1]] %>% filter (rsquared > threshold)
      }
      
      response_names_filt <- perf_filter$response
      
      for (response_value in response_names_filt) {
        
        response_data <- subset(G_target_data_final_df, response == response_value)
        
        target_data_final_df <-  response_data %>% 
          slice_max(order_by = sd_value, n = local_top_var)
        
        plot <- ggplot(target_data_final_df, aes(x = sd_value, y = reorder(var, sd_value))) +
          geom_boxplot() +
          theme_bw() +
          labs(x = "Importance", y = response_value) +
          theme(
            axis.title.x = element_text(size = 8),
            axis.title.y = element_text(size = 8)
          )
        plot_list[[response_value]] <- plot
      }
      
      p2 <- arrangeGrob(grobs = plot_list, plot = FALSE)
      
      combined_plot <- plot_grid(p1, p2, rel_heights = c(1, 0.5), labels = "auto")
      
      #for a taxa of interest
      if(!is.null(taxa)){
        
        # Filter the data for the provided taxa
        target_data <- vi_table %>% dplyr::filter(response == taxa)
        
        # Calculate the mean importance of each variable
        target_data_avg <- target_data %>% 
          dplyr::group_by(var) %>% 
          dplyr::summarise(mean_imp = mean(sd_value)) %>% 
          ungroup()
        
        # Get the top variables
        top_vars <- head(target_data_avg[order(-target_data_avg$mean_imp), ], local_top_var)
        
        # Filter the data to include only the top variables
        target_data_final_df <- target_data %>% dplyr::filter(var %in% top_vars$var)
        
        # Create the boxplot for the target
        p3 <- ggplot(target_data_final_df, aes(x = sd_value, y = reorder(var, sd_value))) +
          geom_boxplot() +
          theme_bw() +
          labs(x = "Importance", y = taxa) +
          scale_y_discrete(label = abbreviate) +
          theme(
            axis.title.x = element_text(size = 8),
            axis.title.y = element_text(size = 8)
          )
        
        #publication ready plot  
        combined_plot <- plot_grid(p1, p2,p3, rel_heights = c(1, 0.5), labels = "auto")
        
      }
      
      
      #for bootstraps
      
      if (!is.null(mrBootstrap_obj)) {
      
      n_response <- ncol(Y)
      complete_df <- cbind(Y, X)
      n_data <- ncol(complete_df)
      
      bind_rows_by_name <- function(list_obj, object_name) {
        filtered_list <- list_obj[names(list_obj) %in% object_name]
        bind_rows(filtered_list)
      }
      
      internal_fit_function <- function(i) {
        
        object_name <- names(complete_df[i])
        
        combined_list <- list()  # Create an empty list to store combined objects
        
        for (j in 1:n_response) {
          combined_object <- map_dfr(mrBootstrap_obj[[j]], bind_rows_by_name, object_name)
          
          if (nrow(combined_object) > 0) {
            
            combined_objfinal <- combined_object #this makes sure it doesnt save the duplicates
            
            combined_list[[j]] <-  combined_objfinal  # Append the combined object to the list
          }
        }
        
        combined_df <- do.call(rbind, combined_list)  # Convert the list to a data frame
        
        return(combined_df)
      }
      
      pd_list <- future_lapply(seq_len(n_data), internal_fit_function, future.seed = TRUE)
      
      
      vi_list <- list()  # Create an empty list to store the plots
      
      for (i in seq_along(pd_list)) {
        
        #extract data for each variable in the models
        df <-  pd_list[[i]]
        
        #calculate the standard deviation of each bootstrap for each taxa
        result <- df %>%
          dplyr::group_by(response, bootstrap) %>%
          summarise(sd_value = sd(value))%>% 
          ungroup() 
        
        #add the variables names back in
        
        df_names <- rep(names(df[1]), nrow(result))
        
        result_update <- cbind(result, var=df_names)
        
        vi_list[[i]] <- result_update
        
      }
      
      vi_df <- do.call(rbind, vi_list)
      #NB#could have a groupCv argument here
      
      G_target_data_avg <- vi_df %>% 
        dplyr::group_by(var) %>% 
        summarise(mean_imp = mean(sd_value)) 
      
      # Create the boxplot for the target
      
      #get top variables
      G_top_vars <- head(G_target_data_avg[order(-G_target_data_avg$mean_imp), ], global_top_var)
      
      G_target_data_final_df <- vi_df %>% 
        dplyr::filter(var %in% G_top_vars$var)
      
      #plot overall importance 
      p1 <- ggplot(G_target_data_final_df, aes(y=reorder(var,sd_value), x=sd_value))+
        geom_boxplot()+
        theme_bw()+
        labs(x="Importance", y="Features")+
        theme(
          axis.title.x = element_text(size = 8),  # Adjust the size as needed
          axis.title.y = element_text(size = 8)
        )
      
      #grid.arrange(p1)
      
      #boxplots for each individual predictor
      
      # Filter the data based on the threshold and calculate mean values for each target
      
      filtered_data <-  vi_df %>%
        dplyr::group_by(response) %>%
        summarise(sd_value = mean(sd_value)) %>%
        right_join(ModelPerf[[1]], by = join_by(response)) %>% 
        filter(roc_AUC > threshold) %>% 
        ungroup()
      
      # Create a list to store the plots
      plot_list <- list()
      
      # Generate boxplots for each target
      for (k in 1:nrow(filtered_data)) {
        
        target <- filtered_data$response[k]###
        
        # Filter the data for the current target
        target_data <-  vi_df %>%
          filter(response == {{target}})
        
        target_data_avg <- target_data %>% 
          dplyr::group_by(var) %>% 
          dplyr::summarise(mean_imp = mean(sd_value)) %>% 
          ungroup()
        
        # Create the boxplot for the target
        
        #get top variables
        top_vars <- head(target_data_avg[order(-target_data_avg$mean_imp), ], local_top_var)
        
        target_data_final_df <- target_data %>% 
          dplyr::filter(var %in% top_vars$var)
        
        plot <- ggplot(target_data_final_df, aes(x = sd_value, y = reorder(var, sd_value))) +
          geom_boxplot() +
          theme_bw()+
          labs(x="Importance", y=target)+
          scale_y_discrete(label=abbreviate)+
          theme(
            axis.title.x = element_text(size = 8),  # Adjust the size as needed
            axis.title.y = element_text(size = 8)
          )
        
        # Add the plot to the list
        plot_list[[k]] <- plot
      }
      
      # Arrange and display the plots in a grid
      #p2 <- grid.arrange(grobs = plot_list, plot = FALSE)
      p2 <- arrangeGrob(grobs = plot_list, plot = FALSE)
      
      #publication ready plot  
      combined_plot <- plot_grid(p1, p2, rel_heights = c(1, 0.5), labels = "auto")
      
      #for a taxa of interest
      if(!is.null(taxa)){

          # Filter the data for the provided taxa
          target_data <- vi_df %>% filter(response == taxa)
          
          # Calculate the mean importance of each variable
          target_data_avg <- target_data %>% 
            dplyr::group_by(var) %>% 
            summarise(mean_imp = mean(sd_value)) %>% 
            ungroup()
          
          # Get the top variables
          top_vars <- head(target_data_avg[order(-target_data_avg$mean_imp), ], local_top_var)
          
          # Filter the data to include only the top variables
          target_data_final_df <- target_data %>% filter(var %in% top_vars$var)
          
          # Create the boxplot for the target
          p3 <- ggplot(target_data_final_df, aes(x = sd_value, y = reorder(var, sd_value))) +
            geom_boxplot() +
            theme_bw() +
            labs(x = "Importance", y = taxa) +
            scale_y_discrete(label = abbreviate) +
            theme(
              axis.title.x = element_text(size = 8),
              axis.title.y = element_text(size = 8)
            )
        
        #publication ready plot  
        combined_plot <- plot_grid(p1, p2,p3, rel_heights = c(1, 0.5), labels = "auto")
        
      }
      
    }
    
    #------------------------------------------------------------------  
    #Importance PCA plot. Responses with similar importance scores group together
    #------------------------------------------------------------------    
 
  vi_table_wide <- vi_table %>% 
    pivot_wider(
                id_cols = 'var',
                names_from = "response",
                values_from = "sd_value") %>% 
                column_to_rownames('var') %>% 
                mutate(across(everything(), ~replace(., is.na(.), mean(., na.rm = TRUE))))#or impute with code below
  
  # #X1 will have missing values. This algorithm will impute them with no impact on the pca
 
       # if(!is.null(X1)) {
  # 
  #   ## First the number of components has to be chosen
  #   ## (for the reconstruction step)
  #   
  # nb <- estim_ncpPCA(vi_table_wide,ncp.max=4)
  #   
  # ## Multiple Imputation
  # vi_table_wide_t <- missMDA::MIPCA( vi_table_wide, ncp = 1)
  # 
  # #get the data
  # vi_table_wide <- vi_table_wide_t$res.imputePCA 
  # 
  # }
  #     
  # else {
  #   
  #   #keep the data as is otherwise
  #   vi_table_wide <- vi_table_wide
  # }
  
      #seems like it doesn't make too much difference
    a.pca <-  t(vi_table_wide) %>% 
      
      prcomp() # do PCA
    
    #-----------------------------------------------------------------------------------------
    #outlier detection
    
    uscores <- a.pca$x %>%
      as.data.frame()
    
    outL <- apply(uscores, 2, function(x) which( (abs(x - median(x)) / mad(x)) > 6 ))
    
    #-----------------------------------------------------------------------------------------
    
    pca_val <-  a.pca  %>%
      tidy(matrix = "eigenvalues")
    
    trans <- t(vi_table_wide)
    
    p3 <- a.pca %>%
      augment(trans) %>% 
      ggplot(aes(.fittedPC1, .fittedPC2)) + 
      geom_point() + 
      ggrepel::geom_label_repel(aes(label = rownames(trans)),
                                box.padding   = 0.35, 
                                point.padding = 0.5,
                                label.size = 0.1,
                                segment.color = 'grey50') +
      theme_bw()
    
    p4 <- a.pca %>%
      tidy(matrix = "eigenvalues") %>%
      ggplot(aes(PC, percent)) +
      geom_col(fill = "#56B4E9", alpha = 0.8) +
      scale_x_continuous(breaks = 1:9) +
      scale_y_continuous(
        labels = scales::percent_format(),
        expand = expansion(mult = c(0, 0.01))
      ) +
      theme_bw() 
    
    combined_plot_PCA <- plot_grid( p3, p4, rel_heights = c(1, 0.5), labels = "auto")
    
    return(list(vi_df,vi_table_wide, combined_plot,  combined_plot_PCA, outLiers=outL,  pca_val=pca_val, scores=uscores))
  }

  