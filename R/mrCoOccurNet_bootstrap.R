#' Generate a MrIML co-occurrence network 
#'
#' This function generates a co-occurrence network from a provided list and calculates strength and directionality of the relationships.
#'
#' @param mrPD_obj A list of model predictions.
#' @param Y The response data.
#'
#' @return A dataframe representing the co-occurrence network with strength and directionality.
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
#'
#'bs_analysis <- mrBootstrap(yhats=yhats_rf,Y=Y, num_bootstrap = 5)
#'assoc_net<- mrCoOccurNet_bootstrap (mrPD_obj=pds , Y=Y)
#'
#'assoc_net_filtered <-  assoc_net %>% 
#'filter(mean_strength > 0.05)

#convert to igraph
#'g <- graph_from_data_frame(assoc_net_filtered, directed=TRUE, vertices=names(Y)) #matching Y data

#'E(g)$Value <- assoc_net_filtered$mean_strength###chnge this as needed
#'E(g)$Color <- ifelse(assoc_net_filtered$direction == "negative", "blue", "red")

# Convert the igraph object to a ggplot object with NMDS layout
#'gg <- ggnetwork(g)

# Plot the graph
#'ggplot(gg, aes(x = x, y = y, xend = xend, yend = yend)) +
#'  geom_edges(aes(color = Color, linewidth = (Value)), curvature = 0.2,
#'             arrow = arrow(length = unit(5, "pt"), type = "closed")) + #makes arrows bigger
#'  geom_nodes(color = "gray", size = degree(g, mode = "out")/2)+#, fill='black', stroke =2) +
#' scale_color_identity() +
#'  theme_void() +
#'  theme(legend.position = "none")  +
#'  geom_nodelabel_repel(aes(label = name),
#'                      box.padding = unit(0.5, "lines"),
#'                       data = gg,
#'                       size=2,
#'                      segment.colour = "black",
#'                      colour = "white", fill = "grey36")} 
#' @export

mrCoOccurNet_bootstrap <- function(mrPD_obj, Y){   #,  variable ='Plas
  
  n_response <- length(Y)
  
  mrPD_obj <- mrPD_obj[[1]] #dont want the plot
  
  PD_yList <- mrPD_obj[seq_along(Y)] #just keep the taxa
  
  result_final_df <- list()
  
  #cooccur_list <- #lapply(PD_yList, function(i){
  for (i in seq_along(PD_yList)) {
    
    #extract data for each variable in the models
    df <-  PD_yList [[i]]
    
    #calculate the standard deviation of each bootstrap for each taxa
    
    df1 <- df %>% mutate(predictor=names(df[1])) %>% 
      rename(class = 1)
    
    result_sd <- df1 %>%
      group_by(target, bootstrap) %>%
      summarize(sd_value = sd(value)) %>% 
      group_by(target) %>% 
      summarize(mean_strength = mean(sd_value), lower_ci = quantile(sd_value, probs = 0.025), upper_ci = quantile(sd_value, probs = 0.975)) %>% 
      mutate(predictor=df1$predictor[1])
    
    
    result_final <- df1 %>%
      spread(class, value) %>%
      mutate(direction = `1` - `0`) %>% 
      group_by(target) %>% 
      summarize(mean_direction=mean(direction)) %>% 
      mutate(direction = ifelse(mean_direction > 0, "positive", "negative")) %>% 
      right_join(result_sd, by='target')
    
    result_final_df[[i]] <- result_final
    
    
  }
  
  complete_net <- do.call(rbind,result_final_df) 
  
  complete_network <- complete_net %>%  rename(taxa_1=predictor,taxa_2= target) %>% 
    dplyr::select(taxa_1, taxa_2, direction, mean_strength, lower_ci, upper_ci, mean_direction)
  
  return(complete_network)
  
}












