
#--------------------------------------------------------------------------------------------------------------
# Use the "ipak" function to install and load multiple R packages. 
# Thanks to Steven Worthington for this code.
# checks to see if packages are installed. Install them if they are not,
# then load them into the R session. Only work for CRAN 
#packages. Details are provided for the others required.

#load packages. Reduce this list.
library(vip)
library(mice)
#library(VIM)
#library(imputeTS)
#library(fastshap)
library(tidymodels)
library(pdp)
library(randomForest)
library(caret)
#library(pROC)
#library(ROCR)
library(missForest)
library(gbm)
#library(iml)
library(tidyverse)
#library(parallel)
library(doParallel)
library(themis)
library(viridis)
library(janitor)
library(hrbrthemes)
library(MRFcov)
library(xgboost)
library(vegan)
library(ggrepel)
library(LEA) 
#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")

#BiocManager::install("LEA")
library(ape)
library(flashlight)
library(devtools)
#install_github("mayer79/flashlight") 
#library(pbapply)
library(flashlight)
# load all function codes. This will disappear when we formally make this a function
source("./R/MrIMLpredicts.R")

source("./R/devianceResids.R")
source(("./R/mrFlashlight.R"))
source("./R/mrProfileplots.R")
#----------------------------------------------------------------------------------------------------------------------



#---------------------------------------------------------------------------------
#REGRESSION DATA FROM FITZPATRICK ET AL Poplar SNP. Proportion of individuals in a population with that SNP
#---------------------------------------------------------------------------------

# read in data file with minor allele freqs & env/space variables
gfData <- read.csv("poplarSNP.ENV.data.4.GF.csv")
envGF <- gfData[,3:13] # get climate & MEM variables
Y <- envGF #for simplicity

# build individual SNP datasets
SNPs_ref <- gfData[,grep("REFERENCE",colnames(gfData))] # reference
GI5 <- gfData[,grep("GI5",colnames(gfData))] # GIGANTEA-5 (GI5)

X <- GI5 #for this example
###############################################################################

#-------------------------------------------------------------------
#Set up the model
#-------------------------------------------------------------------

#user provides the models they would like for each component (model1)

model1 <- #model used to generate yhat
  # specify that the model is a random forest
  linear_reg() %>%
  # select the engine/package that underlies the model
  set_engine("lm") %>%
  # choose either the continuous regression or binary classification mode
  set_mode("regression")

#regression model
model1 <- 
  rand_forest(trees = 100, mode = "regression") %>% 
  set_engine("ranger", importance = c("impurity","impurity_corrected")) %>%
  set_mode("regression")

#for SVM need different tuning paramters. Currently tuning works best for tree-based algorithms. XGR boost wont work on a small dataset like this

#---------------------------------------------------------------------
#Perform the analysis
#---------------------------------------------------------------------  

#models just using features/predictor variables.
yhats <- mrIMLpredicts(X=X,Y=Y, model1=model1, balance_data='no', model='regression', parallel = TRUE) ## in MrTidymodels. Balanced data= up updamples and down downsampled to create a balanced set. For regression there no has to be selected.

#model performance
ModelPerf <- mrIMLperformance(yhats, model1, X=X, model='regression')
ModelPerf[[1]] #predictive performance for individual responses 
ModelPerf[[2]]#overall predictive performance. r2 for regression and MCC for classification

p1 <- as.data.frame(ModelPerf[[1]])

VI <- mrVip(yhats, Y=Y) 
#plot model similarity

#plot variable importance


#not that GLMs in particular wont produce coefficents for features that are strongly colinear and will drop them from the model.
#in this case group cov will have to be changed to reflect features included in the model. 

plot_vi(VI=VI,  X=X,Y=Y, modelPerf=ModelPerf,  cutoff= 0.1, plot.pca='yes', model='regression')#note there are two plots here. PCA is hard to read with > 50 response varianbles





# can build a flaslight object for individual responses 

fl <- mrFlashlight(yhats, X, Y, response = "multi", index=1, model='regression')

plot(light_performance(fl), fill = "orange") +
  
  labs(x = element_blank())

str(fl)

#plot(light_breakdown(fl , new_obs = cbind(X, Y)[1, ]),by = X, v=Y) #prints all responses - need to fix but could be quite handy.

#int <- light_interaction(fl, pairwise=TRUE) #not working, but possible!

#------------------------------------------------------------
#Multiple response

flashlightObj <- mrFlashlight(yhats, X, Y, response = "multi", model='regression')


plot(light_scatter(flashlightObj, v = "bio_1", type = "predicted"))

profileData_pd <- light_profile(flashlightObj, v = "bio_1") #partial dependencies
profileData_ale <- light_profile(flashlightObj, v = "bio_1", type = "ale") #acumulated local effects

#Plot global ALE plot. This the plots the smoothed average ALE value. Have to remove responses that don't respond.

mrProfileplot(profileData_pd , sdthresh =0.01)
mrProfileplot(profileData_ale , sdthresh =0.01)

#------------------------------------------------------------
#Interactions

interactions <-mrInteractions(yhats, X, Y,  model='regression') #this is computationally intensive so multicores are needed. If stopped prematurely - have to reload things

mrPlot_interactions(interactions, X,Y, top_ranking = 20, top_response=5)

save(interactions, 'Fitzpatrick2016interactions')

load('Fitzpatrick2016interactions')
