
#--------------------------------------------------------------------------------------------------------------
# Use the ipak function to install and load multiple R packages. Thanks to Steven Worthington for this code.
# checks to see if packages are installed. Install them if they are not, then load them into the R session. Only work for CRAN 
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
install_github("mayer79/flashlight") 
#library(pbapply)

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

#regression model
model1 <- 
  rand_forest(trees = 100, mode = "regression") %>% 
  set_engine("ranger", importance = c("impurity","impurity_corrected")) %>%
  set_mode("regression")

#for SVM need different tuning paramters. Currently tuning works best for tree-based algorithms.

#---------------------------------------------------------------------
#Perform the analysis
#---------------------------------------------------------------------  

#models just using features/predictor variables.
yhats <- mrIMLpredicts(X=X,Y=Y, model1=model1, balance_data='no', model='regression', parallel = TRUE) ## in MrTidymodels. Balanced data= up updamples and down downsampled to create a balanced set. For regression there no has to be selected.

# can build a flaslight object for individual responses 

fl <- mrFlashlight(yhats, X, Y, response = "single", index=1, mod='regression')

plot(light_performance(fl), fill = "orange") +
  
  labs(x = element_blank())

str(fl)

#plot(light_breakdown(fl , new_obs = cbind(X, Y)[1, ]),by = X, v=Y) #prints all responses - need to fix but could be quite handy.

#int <- light_interaction(fl, pairwise=TRUE) #not working, but possible!

#------------------------------------------------------------
#Multiple response

flashlightObj <- mrFlashlight(yhats, X, Y, response = "multi", mod='regression')

plot(light_profile(flashlightObj, v = "bio_1", type = "ale"))
plot(light_profile(mfl, v = "bio_1", type = "ale"))

#plots
