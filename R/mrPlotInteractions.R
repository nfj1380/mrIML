#'Plots global interactions as well as individual response interaction importance.
#'@param interactions A \code{dataframe} data frame generated from mrInteractioms function 
#'@param Y A \code{dataframe} response data set
#'@param X A \code{dataframe} feature data set
#'@param top_ranking A \code{numeric} determines how many of the strongest feature interacions to view/include
#'@param top_response A \code{numeric} how many of the response variables with the strongest interactions to view
#'@details
#'1st plot: Barplots showing the mean and cumulative importance of each of the top pairs of interactions in the model.
#'2nd plot: Barplot of the responses with the strongest interactions
#'3rd plot: Barplots of the strongest interactions for each of the top response variables.
#'interactions <-mrInteractions(yhats, X, Y) #this is computationally intensive so multicores are needed. If stopped prematurely - have to reload things
#'mrPlot_interactions(Interact, X,Y, top_ranking = 3, top_response=3)
#'@export 


mrPlot_interactions <- function(interactions, X,Y, top_ranking = 3, top_response=10 ){
  
  n_features <- names(X)
  
  variable_interactions <- as.data.frame((t(utils::combn( n_features, m = 2)))) %>% 
    unite('variables', V1:V2, sep='*')
  
  colnames(interactions) <- names(Y)
  
 meanInteractions <- as.data.frame( rowMeans(interactions) ) #calculate average
  names( meanInteractions )[1] <- c('meanInt')
  
  sumInteractions <- as.data.frame( rowSums(interactions) ) #calculate average
  names(sumInteractions )[1] <- c('sumInt')
  
  
  intData <- as.data.frame(interactions)
  intData <- cbind(variable_interactions, intData, meanInteractions, sumInteractions )
  
  inDataOrered <-   intData %>% 
    arrange( variables, meanInt)

  inDataOrered_top <-  inDataOrered[1L:top_ranking, ]
  
    p1 <-  ggplot( inDataOrered_top, aes(y=variables, x= meanInt)) + 
      theme_bw()+
      labs(x= "Mean interaction importance", y='Feature interactions')+
      geom_bar(stat="identity")
    
    p2 <- ggplot(inDataOrered_top, aes(y=variables, x= sumInt)) + 
     theme_bw()+
      labs(x= "Cumulative interaction importance", y='Feature interactions')+
       geom_bar(stat="identity")

      
    grid.arrange(p1,p2, nrow = 1) #plotting both ensures that the cumulative score isn't 
    #biased towards some strong interactions for some predictors
    
    
 #-----------------------------------------------------------------------------------
    
 #select SNPS most effected by interactions for top 10 features
    MostImp <- as.data.frame(colSums(inDataOrered_top[-1]))
    names(MostImp ) <- 'sumInteract' 
    
    MostImp_t <-  MostImp %>% 
      t() %>% 
      as.data.frame()
    
    MostImp_t$meanInt<- NULL #remove these stats as they are not needed
    MostImp_t$sumInt<- NULL
    
    MostImp_f <-  MostImp_t %>% 
      t() %>% 
      as.data.frame() %>% 
      rownames_to_column()
  
      MostImp_ordered <- MostImp_f %>% 
        arrange(desc(sumInteract))
    
      top_int_response <-  as.data.frame(MostImp_ordered[1L:top_ranking, ])
  
      readline(prompt="Press [enter] to continue for response with strongest interactions")
      
  p3 <- ggplot(top_int_response, aes(y=rowname, x= sumInteract)) + #cant get this descending for some reason
        theme_bw()+
        labs(x= "Cumulative interaction importance", y='Response')+
        geom_bar(stat="identity")
  
  print(p3) 
  readline(prompt="Press [enter] to continue for individual response results")   
  t_inDataOrered_top <- as.data.frame(t(inDataOrered_top)) %>% 
        janitor::row_to_names(row_number = 1) %>% 
        rownames_to_column()
      
  topIntC <- filter( t_inDataOrered_top, rowname %in%  top_int_response$rowname)
  
   charvec <- as.data.frame(rep(topIntC$rowname, top_ranking))
    names( charvec) <- 'Response'
  
   topIntC_plotData <- topIntC %>% 
      gather( key ='rowname', value = importance) %>% 
      bind_cols(charvec)
   
   topIntC_plotData$importance <- as.numeric(topIntC_plotData$importance)
   
   
   p4 <- ggplot(topIntC_plotData, aes(fill= rowname , y=importance, x=rowname)) + 
     geom_bar(position="dodge", stat="identity") +
     scale_fill_viridis(discrete = T, option = "E") +
     ggtitle("Individual interaction models") +
     facet_wrap(~Response) +
     theme_ipsum() +
     theme(axis.title.x=element_blank(),
           axis.text.x=element_blank(),
           axis.ticks.x=element_blank())+
     labs(fill='Feature set') 
   
  print(p4)
}
