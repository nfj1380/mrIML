#'Calculates resistance components from a list of pairwise resistance surfaces.
#'@param fileneame A \code{character} this is the location where the resistance surfaces are stored.
#'@param p_val A \code{numeric} this sets the significance threshold for axes in explaining variance in the original resistance matrix based on redundancy analysis. In effect this filters out axes that dont explain variance.
#'@example
# Y <- resist_components(filename = 'Bobcat_cs_matrices', p_val=0.01)
#'@details Outputs a data frame of significant resistance components for each matrix in the target folder. These data can be combined with non-pairwise matrix data.
#'@export 

resist_components <- function (filename = filename, p_val=p_val){
  
  files <- list.files(paste(filename))

  n_matrix <- length(files)
  
#for (i in 1:length(files)){
final_d <- sapply(seq(1, n_matrix), function(i) {

  data_resist <-  read.csv(paste0("./", filename,'/', files[i]))  ### 
  #need to be csv in a folder within your working directory
  data_resist[1] <- NULL #remove row information. It is important that the matrix is symetric
  
  siteData<- as.data.frame(names(data_resist)) #colnames should be names of the population/sites/etc
  
  names(siteData) <- c('Site')#for identifiability

  #names need to match in each matrix
   data_resist <- (data_resist/max(data_resist)) #turns into a dissim matrix
   
   res <- ape::pcoa(data_resist)
   
   l1 <- round(res$values$Relative_eig[1], 2) #variance explained by pcoa 1

   l2 <- round(res$values$Relative_eig[2], 2) #variance explained by pcoa 2 etc
   
   n_axes <- ncol(res$vectors)
   
   pcdat <- as.data.frame(res$vectors)
   
   sigPCs <- NULL
   
#------------------------------------------   

   
    for(j in 1:n_axes){

     vec <- pcdat[j] ####
     
     axisName <- names(vec)
     
     names(vec) <- c('pc_axis')
     
     ## Basic dbRDA Analysis
     vare.cap <- capscale(data_resist ~ vec$pc_axis, parallel = cl) #  sqrt.dist= TRUE avoids negatuve eiganvalues but not working
     
     sig <- anova(vare.cap) #999 permutations
     
     pval <- sig$`Pr(>F)`
     
    # list(axisName = axisName, pval= pval[1])  
     sigPCs[[j]] <- c(axisName, pval[1])  ####
     
   }
   #------------------------------------------
   
   
   PcoA2D_AT <- res$vectors[,1:2]
   # so here simply the Euclidean distance is used, see distances above
   
   PcoA2D_AT <- as.data.frame(PcoA2D_AT)
   
   names(PcoA2D_AT)[1:2] <- c('PCoA1', 'PCoA2')
   
   #plot it
   Tr_PcoA <- ggplot(PcoA2D_AT, aes(x = PCoA1, y = PCoA2, label = siteData$Site))+ #would be good to make this 3D at some point
     geom_point(size =2) +
     geom_text(col = 'black', size=4, check_overlap=TRUE)+ #dont want to print every name
     labs(y = paste("PCoA-2 (", scales::label_percent()(l2), " variance explained)", sep=''), x= paste("PCoA-1 (", scales::label_percent()(l1), " variance explained)", sep=''))+
     theme_bw()+
     ggtitle(paste('PCoA', gsub('.csv','', files[i]))) ###
   print( Tr_PcoA)
   
 # final_d <- NULL
  # all_d <- NULL
   
   sigPC <- as.data.frame(do.call(rbind, sigPCs)) ###
   
   names(sigPC) <-  c('Axis', 'P-value')
   
    sigPC$`P-value` <- as.numeric(as.character(sigPC$`P-value`))
    
   # all_d[[i]] <- sigPC ###
   
   threholdScore <- sigPC  %>%  filter( `P-value`<= p_val)
   
   pcdatT <- as.data.frame(t(pcdat))
  
   ReducedAxis <- pcdatT %>% 
     rownames_to_column() %>% 
     filter (rowname==threholdScore$Axis)
   
   ReducedAxisT <- as.data.frame(t( ReducedAxis)) %>% 
     janitor::row_to_names(row_number = 1) 
    names(ReducedAxisT) <- paste(gsub('.csv','', files[i]), names(ReducedAxisT), sep='_') ###
    
    #ReducedAxisT <-as.numeric(as.character(ReducedAxisT))
    
    ReducedAxisT %>% 
      cbind(siteData)

   
})

final_d <- purrr:: reduce(final_d , left_join, by = 'Site')
row.names(final_d) <- final_d$Site 
final_d$Site <- NULL
final_d <- mutate_all(final_d, function(x) as.numeric(as.character(x)))
final_d


}






