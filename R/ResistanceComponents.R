
#-----------------------------------------------------------
 #sig_pcoa <- function (data_resist, SiteData, spatial){ #GDM version

resist_components <- function (data_resist, SiteData, pval=0.05){
  
   data_resist <- (roadsWS/max(roadsWS)) #turns into a dissim matrix
   
   res <- ape::pcoa(data_resist)
   
   l1 <- round(res$values$Relative_eig[1], 2) #variance explained by pcoa 1
   label_percent()(l1)
   l2 <- round(res$values$Relative_eig[2], 2) #variance explained by pcoa 2 etc
   
   n_axes <- ncol(pcoa$vectors)
   
   pcdat <- as.data.frame(res$vectors)
   
   sigPCs <- NULL
   
#------------------------------------------   
   #sigPC <- sapply(seq(1,n_axes), function(i){
    
     #data_ws<- cbind(SiteData, data_resist )
   
   for(i in 1:n_axes){

     vec <- pcdat[i]
     
     axisName <- names(vec)
     
     names(vec) <- c('pc_axis')
     
     ## Basic dbRDA Analysis
     vare.cap <- capscale(data_resist ~ vec$pc_axis, sqrt.dist= TRUE) #avoids negatuve eiganvalues
     
     sig <- anova(vare.cap)
     
     pval <- sig$`Pr(>F)`
     
    # list(axisName = axisName, pval= pval[1]) 
     sigPCs[[i]] <- c(axisName, pval[1]) 
     
     #ideally this would work with GDM but I get an error with the varImp function for some reason
     
     #pred <- cbind(SiteData,pcdat[i], spatial ) ###
     #gdmTab <- gdm::formatsitepair(data_ws, bioFormat=3,XColumn="Lat", YColumn="Long",
     # siteColumn="Site", predData=pred)
     
     # gdmRes <-  gdm::gdm(gdmTab, geo=F) #dont need geo for now
     
     #mod.test<- gdm::gdm.varImp(gdmTab, geo=FALSE, splines = NULL, knots = NULL, 
     #fullModelOnly = FALSE, nPerm = 100, parallel = TRUE, cores = 2,
     #outFile = NULL)
   }
   #------------------------------------------
   
   sigPC <- as.data.frame(do.call(rbind, sigPCs))

   names(sigPC) = c('Axis', 'P-value')
   
   PcoA2D_AT <- res$vectors[,1:2]
   # so here simply the Euclidean distance is used, see distances above
   
   PcoA2D_AT <- as.data.frame(PcoA2D_AT)
   
   names(PcoA2D_AT)[1:2] <- c('PCoA1', 'PCoA2')
   
   
   Tr_PcoA <- ggplot(PcoA2D_AT, aes(x = PCoA1, y = PCoA2, label = siteData$Site))+
     geom_point(size =2) +
     geom_text(col = 'black', size=4, check_overlap=TRUE)+
     #labs(y= 'PcOA-2 (40% variance explained)', x= 'PcOA-1 (18% variance explained)')+ #these are based on l1/l2
     labs(y = paste("PcOA-2 (", label_percent()(l2), " variance explained)", sep=''), x= paste("PcOA-1 (", label_percent()(l1), " variance explained)", sep=''))+
     theme_bw()
  
   sigPC$`P-value` <- as.numeric(as.character(sigPC$`P-value`))
   
   threholdScore <- test %>%  filter( `P-value`< p_val)
   
   ReducedAxis <- pcdat %>% filter
   
     
   print( Tr_PcoA)
   
   return(sigPC)
   
}
   
