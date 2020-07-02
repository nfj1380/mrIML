#'plot_mrPd: Plots multi-response partial dependencies
# 
#'@param mrpd A \code{list} dataframe from mrPdP.
#'@param Y  A \code{dataframe} #is the feature dataset
#'
#'@details
#'plotting function for mrPdp. 
#'
#'@example
#'plot_mrPd(testPdp, Y)

plot_mrPd <- function (mrpd, Y) {
  
  #create a data frame removing duplicate column names (feature names in this case). Keeps feature values (drop=FALSE)
  m <- as.data.frame(sapply(unique(colnames(testPdp)), function(i) rowMeans(testPdp[,grepl(i, colnames(testPdp)), drop=FALSE])))
  
  n_features <- names(Y)
  
  #create a long dataframe from the wide
  data_l <- gather(m,  key ='Species', value = pd, -c(paste0(n_features)))
  
  p <- ggplot(data_l, aes(x = scale.prop.zos, y = pd)) + #need to work out how to automatically add fetures to x
    geom_line(aes(color = Species), size=1.2) + 
    scale_fill_viridis()+
    theme_bw()

  p
  
  #want functionality to examine a single response variable or the complete set.
  
  #this wont work for models with more than one feature set
}
  
  
  