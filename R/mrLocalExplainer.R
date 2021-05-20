#' Run local explanation methods for individual data points
#'
#' @param dat A \code{dataframe} data frame of predictor variables
#' @param model_dat A \code{workflow} workflow object containing the machine learning model
#' @param outcome \code{vector} vector containing the outcome for each data instance
#'
#' @return 
#' @export
#'
#' @examples
#' #mrLocalExplainer(data, yhats, data$Class)
#' ## will return matrix of phi values, individualized and aggregated plots 
#' 


mrLocalExplainer <- function(dat, model_dat, outcome){

#create list objects needed later  
indiv_plots <- list()
mat1 <- list()

#ensure outcome is stored as a factor
outcome <- as.factor(outcome)

#create colors needed for each outcome
num_levels <- length(levels(outcome)) #number of levels in the outcome factor

my_colors <- pal_uchicago("dark")(9)[1:num_levels] #select same number of colors as there are levels 

#create vector of colors according to outcome level 
color_levels <- as.numeric(outcome) 

for (i in 1:num_levels) {
color_levels[which(color_levels == i)] <- my_colors[i]
}

#define prediction function
predict.function = function(model, new_observation) {
  predict(model, new_observation, "prob")[,2]
}

#pull the model fit from a workflow object
model1 <- pull_workflow_fit(model_dat[[1]]$mod1_k)

#create explainer object
model_explained <- explain.default(model1$fit, dat, predict.function = predict.function)

#determine number of individuals
n <- length(dat[,1])

#breakDown model and individual plots
for (i in 1:n) {
  
  #calculate contribution values 
  f1 = broken(model = model_explained, 
              new_observation = dat[i,], 
              data = dat, 
              predict.function = model_explained$predict_function, 
              keep_distributions = TRUE)

  #create data frame of variables and relative contribution values 
  data1<-data.frame(y = f1$contribution,
                    x = f1$variable) %>%
    filter(!x=="(Intercept)" &!x=="final_prognosis") #remove unnecessary results
  
  #create feature, value and fplot columns
  data1<-data1 %>% mutate(x = paste0(map_chr(str_split(x, "\\+"), 2)))
  data1$feature<-str_split(data1$x,'=', simplify = TRUE)[,1]
  data1$feature<-trimws(data1$feature, which = c("both"))
  data1$svalue<-str_split(data1$x, '=', simplify = TRUE)[,2]
  data1$svalue<-trimws(data1$svalue, which = c("both"))
  data1$feature<-as.character(data1$feature)
  data1$fplot <- paste(data1$feature, "==", data1$svalue)
  
  #create final data frame of contribution results 
  values <- data1$y
  f <- data1$feature
  obs <- data1$svalue
  class <- rep(outcome[i], times = 42)
  ex_data <- cbind(class, f, values, obs)
  mat1[[i]] <- ex_data #store in list object
  
  #produce individual waterfall plots
  indiv_plots[[i]] <- ggplot(data = data1, aes(x = reorder(fplot, -y), y = y)) +
    labs(y = "phi")+
    labs(x = "",subtitle = outcome[i]) +
    geom_bar(stat = "identity", fill = color_levels[i]) + #color will depend on individuals outcome class
    coord_flip() +
    guides(fill = FALSE)+
    theme(text = element_text(size = 14, face = "bold"),
          axis.text.x = element_text(size = 14), 
          axis.text.y = element_text(size = 14, lineheight = 0.7), 
          legend.position = "none") + 
    ggtitle(paste(i)) +
    scale_color_uchicago(palette = "dark") +
    scale_fill_uchicago(palette = "dark")

}

#combine list object into a dataframe
tab1 <- do.call(rbind, mat1[1:n])
tab1 <- as.data.frame(tab1) 
tab1$values <- as.numeric(tab1$values)

tab2 <- tab1 %>% 
  mutate(values_abs = abs(values)) #add absolute values - allows us to observe contributions size regardless of direction

#calculate mean values of variable contributions
tab3 <- aggregate(tab2[,5], list(tab2$f), mean) #aggregates variables and presents the mean values for numeric columns
colnames(tab3)[1] <- "f" 
tab3_ordered <- tab3[order(tab3$x),] #order variables by size
order_only <- tab3_ordered$f #character vector of variable order

#apply variable order to final data frame
tab2$f <- as.factor(tab2$f)
tab2$f <- as.factor(tab2$f, levels = order_only) #not plotting in order at the moment

#summary plot of variable contributions
print(ggplot(tab2, aes(x = f, y = values, fill = as.factor(class))) + 
  geom_boxplot(position = position_dodge(.9)) + 
  geom_hline(aes(yintercept = 0.00), linetype = "dashed") +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16, lineheight = 0.7),
        text = element_text(size = 18, face = "bold"), 
        legend.position = "bottom") + 
  labs(y = "phi", x = NULL) +
  guides(fill=guide_legend(title=outcome)) +
  coord_flip() + 
  scale_color_uchicago() + 
  scale_fill_uchicago())

#return the following objects to the global environment
LE_matrix <<- mat1
LE_indiv_plots <<- indiv_plots 

}