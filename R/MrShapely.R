#' Generate SHAP (SHapley Additive exPlanations) Plots for Multiple Models and Responses
#'
#' This function generates SHAP (SHapley Additive exPlanations) plots for multiple models and responses.
#'
#' @param yhats A list of model prediction objects. Each object should contain a model, data, and class information.
#' @param MultRespVars A data frame containing response variables for the prediction.
#' @param taxa An optional vector specifying which responses to include based on their indices.
#' @param x_features A character vector specifying the features to consider in the plots.
#' @param y_features A character vector specifying the response features for interaction plots.
#' @param kind A character string specifying the type of plots (e.g., "beeswarm" for feature effect plot, "bar" for variable importance plot, or "both").
#' @param max_display An integer specifying the maximum number of features to display.
#' @param interactions A logical indicating whether to create interaction effect plots.
#' @param color_var A variable to use for coloring in the dependency plots.
#' @param getFeaturePlot A logical indicating whether to generate feature effect plots.
#' @param getDependencyPlot A logical indicating whether to generate dependency plots.
#' @param getInteractionPlot A logical indicating whether to generate interaction plots.
#' @param num_cores An integer specifying the number of CPU cores to use for parallel processing.
#' @param class_selection An optional vector specifying which classes to include in the plots.
#'
#' @return A ggplot object containing SHAP plots for the specified responses and features. Note that this function may not work for some algorithm classe (e.g., neural nets)
#'
#' @examples
#' # Example usage:
#' MrShapely(yhats, MultRespVars = Resp)
#' @export
#' 

MrShapely <- function(yhats, MultRespVars = Resp,
                      taxa = NULL,
                      kind = "beeswarm",
                      max_display = 15L,
                      color_var = NULL, getFeaturePlot = TRUE,
                      getDependencyPlot = TRUE, get2DDependencyPlot = TRUE,
                      num_cores = 2,
                      class_selection = NULL) {
  
  # Parallelization
  plan(future::multisession, workers = num_cores)
  
  # Extract model, data, and response information
  mod_list <- future_lapply(yhats, function(yhat) yhat$mod1_k,future.seed = TRUE)
  model_list <-lapply(mod_list, extract_fit_parsnip)#This part does not seem to work with future_lapply
  Xdata_list <- future_lapply(yhats, function(yhat) as.data.frame(yhat$data),future.seed = TRUE)
  X_i_list <- future_lapply(yhats, function(yhat) as.data.frame(select(yhat$data, -class)),future.seed = TRUE)
  Y_i_list <- future_lapply(yhats, function(yhat) yhat$data$class,future.seed = TRUE)
  
  # Calculate Shap_kernel and shapobj for all responses using future_lapply
  shapobj_list <- future_lapply(1:length(model_list), function(i) {
    model <- model_list[[i]]
    X_i <- X_i_list[[i]]
    Xdata <- Xdata_list[[i]]
    Y_i <- Y_i_list[[i]]
    
    # Deals with classification cases
    if (model$spec$mode == "classification") {
      predfun <- function(model, newdata) {
        preds <- predict(model, as.data.frame(newdata), type = 'response') # probabilities
        return(cbind(1 - preds, preds)) # for both classes
      }
      
      if (inherits(model$fit, "glm")) {
        Shap_kernel <- kernelshap(model$fit, X = X_i,
                                  pred_fun = predfun,
                                  bg_X = Xdata)
        shapobj <- shapviz(Shap_kernel)
        names(shapobj) <- levels(Y_i)
      } else if (inherits(model$fit, "randomForest") || inherits(model$fit, "ranger")) {
        Shap_kernel <- kernelshap(model, X = X_i,
                                  bg_X = Xdata,
                                  type = "prob")
        shapobj <- shapviz(Shap_kernel)
        names(shapobj) <- levels(Y_i)
      } else {
        shapobj <- shapviz(model$fit,
                           X_pred = data.matrix(X_i),
                           X = Xdata)
        names(shapobj) <- levels(Y_i)
      }
      
      return(shapobj)
    }
    
    # Deals with regression cases
    if (model$spec$mode == "regression") {
      if (inherits(model$fit, "lm")) {
        Shap_kernel <- kernelshap(model$fit, X = X_i,
                                  bg_X = Xdata)
        shapobj <- shapviz(Shap_kernel)
      } else if (inherits(model$fit, "randomForest") || inherits(model$fit, "ranger")) {
        Shap_kernel <- kernelshap(model, X = X_i,
                                  bg_X = Xdata)
        shapobj <- shapviz(Shap_kernel)
      } else {
        shapobj <- shapviz(model$fit,
                           X_pred = data.matrix(X_i))
      }
      
      return(shapobj)
    }
  },future.seed = TRUE)
  
  
  # Subset the results based on the 'taxa' parameter
  if (!is.null(taxa)) {
    if (is.numeric(taxa)) {
      if (length(taxa) == 1) {
        # print("Subset shapobj_list based on 'taxa'.")
        shapobj_list <- list(shapobj_list[[taxa]])
        ResponseNames <- colnames(MultRespVars)[taxa]
      } else if (length(taxa) > 0 && all(taxa > 0 & taxa <= length(shapobj_list))) {
        # print("Subset shapobj_list based on 'taxa'.")
        shapobj_list <- shapobj_list[taxa]
        ResponseNames <- colnames(MultRespVars)[taxa]
      } else {
        stop("Invalid value for 'taxa'. Please provide valid indices.")
      }
    } else {
      stop("Invalid value for 'taxa'. Please provide numeric indices.")
    }
  } else {
    ResponseNames <- colnames(MultRespVars)
    taxa <- seq_along(shapobj_list)
  }
  
  # Initialize an empty list to store plots with labels
  plots_with_labels <- list()
  
  # [1] FEATURE EFFECT PLOT FUNCTION
  if (isTRUE(getFeaturePlot)) {
    if (is.null(taxa)) {
      taxa_to_iterate <- seq_along(shapobj_list)
    } else {
      taxa_to_iterate <- seq_along(shapobj_list)
    }
    s
    feature_plots_with_labels <- future_lapply(taxa_to_iterate, function(i) {
      response_name <- ifelse(is.null(ResponseNames[i]), "", ResponseNames[i])
      
      if (model_list[[i]]$spec$mode == "regression") {
        label <- response_name
        shapobj <- shapobj_list[[i]]
        plot_obj <- sv_importance(shapobj, kind = kind, max_display = max_display) + labs(title = label)
      } else {
        class_list <- shapobj_list[[i]]
        print(paste("Index:", i, "Length of shapobj_list:", length(shapobj_list)))
        class_feature_plots <- future_lapply(names(class_list), function(class_name) {
          if (!is.null(class_selection) && !(class_name %in% class_selection)) {
            return(NULL)
          }
          
          class_obj <- class_list[[class_name]]
          label <- ifelse(is.null(class_name), response_name, paste(response_name, class_name, sep = " - "))
          plot_obj <- sv_importance(class_obj, kind = kind, max_display = max_display) + labs(title = label)
          return(plot_obj)
        },future.seed = TRUE)
        
        class_feature_plots <- class_feature_plots[!future_sapply(class_feature_plots, is.null)]
        
        if (length(class_feature_plots) > 0) {
          return(class_feature_plots)
        } else {
          return(NULL)
        }
      }
    },future.seed = TRUE)
    
    plots_with_labels <- c(plots_with_labels, feature_plots_with_labels)
  }
  # END OF FEATURE EFFECT PLOT FUNCTION
  
  # [2] DEPENDENCY PLOT FUNCTION
  if (isTRUE(getDependencyPlot)) {
    if (is.null(taxa)) {
      taxa_to_iterate <- seq_along(shapobj_list)
    } else {
      taxa_to_iterate <- seq_along(shapobj_list)
    }
    
    dependency_plots_with_labels <- 
      future_lapply(taxa_to_iterate, function(i) {
        response_name <- ifelse(is.null(ResponseNames[i]), "", ResponseNames[i])
        
        if (model_list[[i]]$spec$mode == "regression") {
          label <- response_name
          shapobj <- shapobj_list[[i]]
          colnames_shapobj <- colnames(shapobj)
          dependency_plots <- future_lapply(colnames_shapobj, function(x_feature) {
            plot_obj <- sv_dependence(shapobj, x_feature, color_var = NULL) + labs(title = label)
            return(plot_obj)
          })
          ncol <- length(colnames_shapobj)
          arranged_plots <- do.call(gridExtra::grid.arrange, c(dependency_plots, ncol = ncol))
          return(arranged_plots)
        } else {
          class_list <- shapobj_list[[i]]
          class_dependency_plots <- future_lapply(names(class_list), function(class_name) {
            if (!is.null(class_selection) && !(class_name %in% class_selection)) {
              return(NULL)
            }
            class_obj <- class_list[[class_name]]
            label <- ifelse(is.null(class_name), response_name, paste(response_name, class_name, sep = " - "))
            colnames_class_obj <- colnames(class_obj)
            dependency_plots <- future_lapply(colnames_class_obj, function(x_feature) {
              plot_obj <- sv_dependence(class_obj, x_feature, color_var = NULL) + labs(title = label)
              return(plot_obj)
            }, future.seed = TRUE)
            ncol <- length(colnames_class_obj)
            arranged_plots <- do.call(gridExtra::grid.arrange, c(dependency_plots, ncol = ncol))
            return(arranged_plots)
          }, future.seed = TRUE)
          
          class_dependency_plots <- class_dependency_plots[!future_sapply(class_dependency_plots, is.null)]
          
          if (length(class_dependency_plots) > 0) {
            return(class_dependency_plots)
          } else {
            return(NULL)
          }
        }
      }, future.seed = TRUE)
    
    plots_with_labels <- c(plots_with_labels, dependency_plots_with_labels)
  }
  
  # [3] INTERACTION EFFECT PLOT FUNCTION
  if (isTRUE(get2DDependencyPlot)) {
    if (is.null(taxa)) {
      taxa_to_iterate <- seq_along(shapobj_list)
    } else {
      taxa_to_iterate <- seq_along(shapobj_list)
    }
    
    
    interaction_plots_with_labels <- 
      future_lapply(taxa_to_iterate, function(i) {
        response_name <- ifelse(is.null(ResponseNames[i]), "", ResponseNames[i])
        
        if (model_list[[i]]$spec$mode == "regression") {
          label <- response_name
          shapobj <- shapobj_list[[i]]
          colnames_shapobj <- colnames(shapobj)
          interactions <- combn(colnames_shapobj, 2, simplify = FALSE) # Get all combinations of features
          interaction_plots <- lapply(interactions, function(int) {
            x_feature <- int[1]
            y_feature <- int[2]
            interaction_label <- paste(response_name, sep = " - ")
            interaction_plot <- sv_dependence2D(shapobj, x_feature, 
                                                y_feature) + labs(title = interaction_label)
            return(interaction_plot)
          })
          return(interaction_plots)
        } else {
          class_list <- shapobj_list[[i]]
          interaction_plots <- list()
          
          for (class_name in names(class_list)) {
            if (!is.null(class_selection) && !(class_name %in% class_selection)) {
              next
            }
            class_obj <- class_list[[class_name]]
            colnames_class_obj <- colnames(class_obj)
            interactions <- combn(colnames_class_obj, 2, simplify = FALSE) # Get all combinations of features
            interaction_plots_for_feature <- lapply(interactions, function(int) {
              x_feature <- int[1]
              y_feature <- int[2]
              interaction_label <- paste(response_name, class_name, sep = " - ")
              interaction_plot <- sv_dependence2D(class_obj, x_feature, y_feature) + labs(title = interaction_label)
              return(interaction_plot)
            })
            interaction_plots <- c(interaction_plots, interaction_plots_for_feature)
          }
          
          if (length(interaction_plots) > 0) {
            return(interaction_plots)
          } else {
            return(NULL)
          }
        }
      })
    
    plots_with_labels <- c(plots_with_labels, interaction_plots_with_labels)
  }
  
  # Flatten the list
  plots_with_labels <- unlist(plots_with_labels, recursive = FALSE)
  
  # Remove NULL elements
  plots_with_labels <- plots_with_labels[!sapply(plots_with_labels, is.null)]
  
  # Stop parallel processing
  future::plan(future::sequential)
  
  ncol <- max(1, ceiling(sqrt(length(plots_with_labels))))
  
  # Create the final plot by stacking responses
  final_plot <- do.call(gridExtra::grid.arrange, c(plots_with_labels, ncol = ncol))
  
  final.plot <- as.ggplot(final_plot)
  
  
  if (getFeaturePlot) {
    if (kind == "beeswarm") {
      final.plot <- final.plot +
        labs(title="Feature Effect Plot")+
        theme(plot.title.position = "plot")
    } else if (kind == "bar") {
      final.plot <- final.plot +
        labs(title="Feature Importance Plot")+
        theme(plot.title.position = "plot")
    } else if (kind == "both") {
      final.plot <- final.plot +
        labs(title="Feature Effect and Importance Plot")+
        theme(plot.title.position = "plot")
    }
  }
  
  if (getDependencyPlot) {
    final.plot <- final.plot +
      labs(title="Dependency Plot")+
      theme(plot.title.position = "plot")
  }
  
  if (get2DDependencyPlot) {
    final.plot <- final.plot +
      labs(title="Interaction Effect Plot")+
      theme(plot.title.position = "plot")
  }
  
  return(final.plot)
}