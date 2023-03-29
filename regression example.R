library(mrIML)

#other package we need:
library(vip); library(tidymodels);
library(randomForest);library(gbm);
library(tidyverse);library(parallel); 
library(doParallel); library(themis);
library(viridis); library(janitor);
library(hrbrthemes);  library(vegan);
library(flashlight);library(iml);
library(ggrepel); library(ranger);
library(future.apply); library(plyr)

cl <- parallel::makeCluster(3)
future::plan(cluster, workers=cl)

X<-as.data.frame(gfdf[,2:13])
Y<-gfdf[,seq(14,61822,100)]
model_rf <- rand_forest(trees = 10, mode = "regression", mtry = tune(), min_n = tune()) %>% set_engine("randomForest")
yhats_rf <- mrIMLpredicts(X=X,Y=Y, Model=model_rf,
                          balance_data='no', mode='regression', 
                          tune_grid_size=5, seed = sample.int(1e8, 1)
)
ModelPerf <- mrIMLperformance(yhats_rf, Model=model_rf, Y, mode='regression')


                          