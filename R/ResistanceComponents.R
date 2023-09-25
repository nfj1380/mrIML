#'Calculates resistance components from a list of pairwise resistance surfaces.
#'@param foldername A \code{character} this is the location where the resistance surfaces are stored.
#'@param p_val A \code{numeric} this sets the significance threshold for axes in explaining variance in the original resistance matrix based on redundancy analysis. In effect this filters out axes that dont explain variance.
#'@example
# Y <- resist_components(filename = 'Bobcat_cs_matrices', p_val=0.01)
#'@details Outputs a data frame of significant resistance components for each matrix in the target folder. These data can be combined with non-pairwise matrix data.
#'@export 

resist_components <- function (foldername = foldername, p_val=p_val){
  
  files <- list.files(paste(foldername))

  n_matrix <- length(files)
  
#for (i in 1:length(files)){
final_d <- lapply(seq(1, n_matrix), function(i) {

  data_resist <-  read.csv(paste0("./", foldername,'/', files[i]))  ### 
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
   
   
##########Plot##################
   
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
   
  
#------------------------------------------ 
#Select only significant axes
#------------------------------------------   
   
   sigPCs <- NULL
   
    for(j in 1:n_axes){

       vec <- pcdat[j] ####
     
      axisName <- names(vec)
     
       names(vec) <- c('pc_axis')
     
     ## Basic dbRDA Analysis
      vare.cap <- vegan::capscale(data_resist ~ vec$pc_axis, parallel = cl) #  sqrt.dist= TRUE avoids negatuve eiganvalues but not working
     
       sig <- anova(vare.cap) #999 permutations
     
       pval <- sig$`Pr(>F)`
     
     sigPCs[[j]] <- c(axisName, pval[1])  ####
     
   }

   sigPC <- as.data.frame(do.call(rbind, sigPCs)) ###
   
   names(sigPC) <-  c('Axis', 'P-value')
   
    sigPC$`P-value` <- as.numeric(as.character(sigPC$`P-value`))
   
    #select significant pcs
   threholdScore <- sigPC  %>%  filter( `P-value`<= p_val)
   
#------------------------------------------   
   
  #select informative axes
    ReducedAxis <- pcdat %>% ##
      subset(select=threholdScore$Axis)
   
    # add site data
    ReducedAxisDF <-  cbind(siteData,  ReducedAxis)
    rownames(ReducedAxisDF) <- ReducedAxisDF$Site
    ReducedAxisDF$Site <-NULL
    
    #add resistance surface names
    names(ReducedAxisDF) <- paste(gsub('.csv','', files[i]),
                                 names(ReducedAxisDF), sep='_')
    ReducedAxisDF
})
#save as a data frame
   final_d <- as.data.frame(do.call(cbind, final_d))


}






