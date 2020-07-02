#This is an old copy of our previous work - elements maybe useful. DON'T work on this code further. 

#' Perform mutlivariate ensemble learning with models available in caret
#'
#' This function allows you to perform ensemble learning on mutliple response variables
#' @param X is a data.frame with rows as sites or individuals/populations and columns as loci or species or OTUs.
#' @param Y is a data.frame with rows as sites or individuals/populations and columns as features e.g. environmental or host or demographic variables.
#' @param grid is the model specific grid of values to tune each model e.g. random forest. Not quite working properly yet
#' @keywords multivariate ensemble learning
#' @export

mcaret <- function(X, Y, grid ){
  
  if (!is.data.frame(Y)) 
    Y <- data.frame(Y)
  if (!is.data.frame(X)) 
    X <- as.data.frame()
  
  
  n <- nrow(X)
  response.names <- names(X)
  site.names <- colnames(X)
  
  model_results <- NULL
  modIml_list  <- NULL
  
  #performance_results <- data.frame(matrix(vector(), 40,4, 
  #dimnames=list(c(), c("AUC", "Specificity", "Sensitivity", "MCC"))))
  performance_results <- c()
  
  for(i in 1:length(X)) {
    
    # ## create folds for CV
    data <- cbind(X[i], Y)
    colnames(data)[1]=c('class')
    
    
    inTrain <- createDataPartition(y=data$class, p = .8, list = FALSE) #.8 of the dataset for training is reccomended 
    
    data.train.descr <- data[inTrain,-1]
    data.train.class <- data[inTrain,1]
    data.test.descr <- data[-inTrain,-1]
    data.test.class <- data[-inTrain,1]
    
    myFolds <- createMultiFolds(y=data.train.class,k=10,times=10)
    
    #make sure can do in parallel - leads to bugs
    cluster <- makeCluster(detectCores() -1) #leave one for the OS
    registerDoParallel(cluster)
    
    myControl <- trainControl(## 10-fold CV)
      method = "repeatedcv",
      number = 10,
      repeats = 10,
      index = myFolds,
      savePredictions=TRUE,        
      classProbs = TRUE,
      summaryFunction = twoClassSummary,  
      allowParallel=TRUE,
      selectionFunction = "best")
    
    model<- train(data.train.descr,data.train.class, #this bit needs fixing
                  method = "rf",#print(modeltype), #this is where we select which model to use for prediction (out of the 237 available)
                  metric = "ROC",# Reciever opperating characteristic used to test model performance     (sensitivity/specificity)
                  verbose = FALSE,
                  trControl = myControl,
                  tuneGrid = grid, # Optimized parameter tuning is done usig the 'rf.grid' set up
                  allowParallel=TRUE,
                  verboseIter=TRUE)
    
    #stop parallel computing
    stopCluster(cluster)
    registerDoSEQ()
    
    model_results <-c(model_results, list(model))
    sapply(model_results, class)
    
    Y1 <-data[-which(names(data) == "class")] #load data again for the visualization. Maybe remove here.
    X1 <- as.data.frame(data$class)
    
    modIml <-Predictor$new(model, data = Y1, y = X1) #this needs to be fixed
    
    modIml_list <- c(modIml_list , list(modIml))
    
    
    #model performance
    oof <- model$pred 
    oof <- oof[oof$mtry==model$bestTune[,'mtry'],]# consider only the best tuned model to use. Is this is the same for model trypes
    
    repeats<-model$control$repeats
    
    model.performance.cv <- EstimatePerformanceCV(oof=oof,repeats=repeats) #Matthew's correlation coefficient (MCC)
    performance_results <- rbind(model.performance.cv$mean.metrics, performance_results)
  }
  
  names(model_results) <- response.names
  names(modIml_list) <- response.names
  row.names(performance_results) <- site.names
  #plots
  perf <- as.data.frame(performance_results)
  perf <- rownames_to_column(perf, "Response")
  
  perf <- arrange(perf, desc(MCC))
  perf$Response <- factor(perf$Response, levels=perf$Response  )
  Global_avg <- perf %>% summarize(meanMCC = mean(MCC), meanAUC=mean(ACC), meanSEN= mean(SEN), meanSpe = mean(SPE))
  
  pMCC <- ggplot(perf, aes(y=MCC, x=Response))+
    geom_point(stat='identity')+theme_bw()+theme(axis.text = element_text(angle=90,hjust=1))+ggtitle('modeltype')
  print(pMCC)
  
  return(list(model_results, performance_results,Global_avg, modIml_list))
}
