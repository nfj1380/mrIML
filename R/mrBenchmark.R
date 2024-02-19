#' Compare and benchmark disease outbreak risk among and within groups
#'
#' @param data A \code{character} object name of data frame 
#' @param Y A \code{character} column name of variable containing outcome 
#' @param pred A \code{character} column name of variable containing model predicted values
#' @param group A \code{character} column name of variable that individuals should be grouped by
#' @param type A \code{character} specify within group "internal" or among group "external" benchmarking 
#' @param label_by A \code{character} column name of variable representing the individual units. If stated, these will be labeled on the ggplot. By default labels will not be included
#'@examples 
#'mB <- mrBenchmark(data=data, Y='class')
#' @export


mrBenchmark <- function(data = "data", Y = "class", pred = "predicted", group = "group1", label_by = "ID", type = "internal"){
  
  #create objects from character strings
  data1 <- as.data.frame(eval(parse(text=data)))
  outcome1 <- as.factor(eval(parse(text=paste("data1$", Y, sep=""))))
  pred1 <- eval(parse(text=paste("data1$", pred, sep="")))
  group1 <- as.factor(eval(parse(text=paste("data1$", group, sep=""))))
  
  c <- discretize(as.numeric(pred1), cuts = 3, labels = c("Low risk", "Medium risk", "High risk"), keep_na = FALSE, infs = FALSE)
  
  ##among group benchmarking
  if(type == "external"){
    
    print(ggplot(data = data1, aes(x="", y=pred1, color=as.factor(outcome1))) +
            geom_boxplot(outlier.size=-1, lwd = 1.3)+
            facet_grid(as.formula(paste(".~", group)))+ #facetting variable
            scale_y_continuous(limits = c(0, 1))+
            ylab(pred)+
            geom_hline(aes(lty="High risk",yintercept= c$breaks[3])) +
            geom_hline(aes(lty="Low risk",yintercept= c$breaks[2])) +
            scale_linetype_manual(name="Risk level",values=c(1,2)) +
            theme(plot.title = element_text(hjust = 0.5)) +
            theme(text = element_text(size = 17, face = "bold"),
                  axis.title.x=element_blank(),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank()) +
            scale_color_uchicago(palette = "dark", name = Y))
    
  }
  
  if(type == "internal"){
    
    if(is.null(label_by) == FALSE){ 
      
      label1 <- as.factor(eval(parse(text=paste("data1$", label_by, sep=""))))
      
      #determine number of groups
      n <- levels(group1)
      
      #plot per group
      for(i in n){
        
        #only consider data from specified group
        data2 <- data1%>%
          filter(group1 == i)
        
        #create more objects
        outcome2 <- as.factor(eval(parse(text=paste("data2$", Y, sep=""))))
        pred2 <- eval(parse(text=paste("data2$", pred, sep="")))
        label2 <- as.factor(eval(parse(text=paste("data2$", label_by, sep=""))))
        
        #internal plots
        print(ggplot(data2, aes(x=outcome2, y=pred2, color = as.factor(outcome2))) +
                geom_point(aes(color = as.factor(outcome2))) +
                xlab(Y) +
                ylab(pred) +
                geom_hline(aes(lty="High risk",yintercept= c$breaks[3])) +
                geom_hline(aes(lty="Low risk",yintercept= c$breaks[2])) +
                scale_linetype_manual(name="Risk level",values=c(1,2)) +
                ggtitle(i) +
                guides(colour = guide_legend(title = Y)) +
                geom_label_repel(aes(label=label2),
                                 fontface = "bold",
                                 force =10,
                                 box.padding = unit(0.35, "lines"),
                                 point.padding = unit(0.5, "lines"),
                                 arrow = arrow(length = unit(0.01, "npc"), 
                                               type = "open", ends = "last"),
                                 size = 5, max.overlaps = 20, show.legend = F) +
                theme(text = element_text(size = 12, face = "bold"),
                      plot.title = element_text(size = 18, hjust = 0.5),
                      axis.line = element_line(size=0.8, colour = "black"), 
                      axis.text.x = element_text(colour="black", size = 18), 
                      axis.text.y=element_text(colour="black", size = 18), 
                      axis.title.x = element_text(size = 18), 
                      axis.title.y = element_text(size = 18),
                      strip.background = element_rect(color="black",size=1.5, linetype="solid"))+
                scale_color_uchicago(palette = "dark") +
                scale_fill_uchicago(palette = "dark"))}
    }
    
    if(is.null(label_by) == TRUE){
      
      #determine number of groups
      n <- levels(group1)
      
      #plot per group
      for(i in n){
        
        #only consider data from specified group
        data2 <- data1%>%
          filter(group1 == i)
        
        #create more objects
        outcome2 <- as.factor(eval(parse(text=paste("data2$", Y, sep=""))))
        pred2 <- eval(parse(text=paste("data2$", pred, sep="")))
        
        #internal plots
        print(ggplot(data2, aes(x=outcome2, y=pred2, color = as.factor(outcome2))) +
                geom_point(aes(color = as.factor(outcome2))) +
                xlab(Y) +
                ylab(pred) +
                ggtitle(i) +
                geom_hline(aes(lty="High risk", yintercept = c$breaks[3])) +
                geom_hline(aes(lty="Low risk", yintercept = c$breaks[2])) +
                scale_linetype_manual(name="Risk level",values=c(1,2)) +
                guides(colour = guide_legend(title = Y)) +
                theme(text = element_text(size = 12, face = "bold"),
                      plot.title = element_text(size = 18, hjust = 0.5),
                      axis.line = element_line(size=0.8, colour = "black"), 
                      axis.text.x = element_text(colour="black", size = 18), 
                      axis.text.y=element_text(colour="black", size = 18), 
                      axis.title.x = element_text(size = 18), 
                      axis.title.y = element_text(size = 18),
                      strip.background = element_rect(color="black",size=1.5, linetype="solid"))+
                scale_color_uchicago(palette = "dark") +
                scale_fill_uchicago(palette = "dark"))
      }
      
    }
  }
}