#'Testing functions out on real data.I will turn this into a vignette
#' to go with the paper. Currently we have two test datasets. 
#'One coinfection data set where columns represent parasite species.
#' Rows reflect sites. This dataset has only one feature (scale prop.
#' Currently just set up for binary classification. But multiclass/regression 
#' would be useful additions and easy
#'in the tidymodel environment. Multiclass could be pick SNP identities which would be useful.
#'
#
#lean toward mrIML as a name for the package.

#I'm thinking for the paper having three test datasets plus simulations.
#1. FIV bobcat data (host/pathogen example) 
#2. Amphibian SNP data (landscape genetics).
#3. Moose microbiome data. This would hit the key audience for Molecular Ecology resources.
#4.PRRS data on multi-strain circulation in swine https://www.biorxiv.org/content/10.1101/2020.04.09.034181v2

#There then could be separate papers on confection data extending 
#the approach to focus on co-occurence  (Lion/avian malaria data) 
#as well as a paper tracking between farm transisiion (Gustavo's data)

#load packages. Reduce this list.
library(vip)
library(mice)
library(VIM)
library(imputeTS)
library(fastshap)
library(tidymodels)
library(pdp)
library(randomForest)
library(caret)
library(pROC)
library(ROCR)
library(plyr)
library(missForest)
library(tidyverse)
library(gbm)
library(iml)
library(tidyverse)
library(parallel)
library(doParallel)
library(themis)
library(viridis)
library(janitor)
library(hrbrthemes)
library(MRFcov)
library(xgboost)

# load all function codes. This will disappear when we formally make this a function
source("./R/filterRareCommon.R")
source("./R/MrTidyModels.R")
source("./R/StackPredictions.R")
source("./R/devianceResids.R")
source("./R/filterRareCommon.R")
source("./R/MrtTidyPerf.R")
source("./R/mrvip.R")
source("./R/plot_vi.R")


#Nick C - this and the function below are the functions with tidy model code I made from your original
source("./R/stacked_preds.R") #this does the stacking
source("./R/response_covariance") #this should create the covatiance matrix

#---------------------------------------------------------------------------------

#Viral SNP test data
set.seed(123)
#load data for training
Responsedata <-read.csv("GF SNPs posBinary.csv")
rownames(Responsedata) <-Responsedata$X
Responsedata$X <- NULL
Responsedata[Responsedata==4] <-0 #fix random 4s in response data. 

#feature set
Features <-read.csv("Landscape and host data.csv", row.names = 1, head=T)
#summary(Features)

# instead we can devlop a cut % for interpolation in the func_load_data. Seems to be a random row somehow added to the feature data.
#Features<-na.interpolation(Features, option = "spline")

# # or remove NAs from the feature/predictor data.
FeaturesnoNA<-Features[complete.cases(Features), ];str(Features) #dropping NAs

Y <- FeaturesnoNA #for simplicity

#---------------------------------------------------------------------------------
#coinfection test data # from MRFcov. This data can also be used in this pipeline replacing X/Y
X <- select(Bird.parasites, -scale.prop.zos) #response variables eg. SNPs, pathogens, species....
Y <- select(Bird.parasites, scale.prop.zos) # feature set
#---------------------------------------------------------------------------------

#Optional: Filter rare/common SNPs or species

fData <- filterRareCommon (Responsedata, lower=0.35, higher=0.25) #this removes all SNPs 
X <- fData #for simplicity when comparing
 
#that occur in less than 35% of individuals and > 75% of individuals
#fData <- rownames_to_column(fData, "Individual")#get individual id back
#user can specify any tidy model here. 

#---------------------------------------------------------------------
#Set up the model
#-------------------------------------------------------------------

#user provides the models they would like for each component (model1)
 
model1 <- #model used to generate yhat
  # specify that the model is a random forest
  logistic_reg() %>%
  # select the engine/package that underlies the model
  set_engine("glm") %>%
  # choose either the continuous regression or binary classification mode
  set_mode("classification")

#random forest. Note that these models are tuned   - need to fix this.
 
model1 <- 
  rand_forest(trees = 100, mode = "classification") %>%
  set_engine("ranger", importance = "impurity") %>%
  set_mode("classification")

#Boosted model (xgrboost). Not tuned either currently. This doesn't work with the small bobcat dataset. Much too small
model1 <- 
boost_tree() %>%
  set_engine("xgboost") %>%
  set_mode("classification")

#---------------------------------------------------------------------
#Perform the analysis
#---------------------------------------------------------------------  
  
#models just using features/predictor variables.
yhats <- mrIMLpredicts(X=X,Y=Y, model1=model1, balance_data='no') #model 1 has to be the model used in mrIMLpredicts. Im sure we could fix this

#ideally we shopuld get this to work on the stacked or not stacked models
ModelPerf <- mrIMLperformance(yhats, model1, X=X) #ROC wont work for some reason. But MCC is useful (higher numbers = better fit)
#doesn't work with upsampled data - something to do with mcc? Still unclear
ModelPerf #nice to have ROC curves for each response (each response with a differing colour. Add prevalence

save(ModelPerf, "randF.RData") #save to compare later.

#we can look at variable importance 
VI <- mrVip(yhats, Y=Y)

#plot variable importance

#for interpretation group features. Annoying but I cant see an easier way. 
groupCov <- c(rep ("Host_characteristics", 1),rep("Urbanisation", 3), rep("Vegetation", 2), rep("Urbanisation",1), rep("Spatial", 2), 
              rep('Host_relatedness', 6),rep ("Host_characteristics", 1),rep("Vegetation", 2), rep("Urbanisation",1))  

plot_vi(VI=VI,  X=fData,Y=FeaturesnoNA, modelPerf=ModelPerf, groupCov, cutoff= 0.5) #note there are two plots here.#I get a strange error running this some times 'prop not found'. 
#First plot is overall importance and the second is individual SNP models.
#warning are about x axis labels so can ignore 

testPdp <- mrPdP(yhats, X=X,Y=Y, Feature='Longitude')
#something to do with the expand grid function within partial. Making a vector too long?
#invalid 'times' value 

#indiv SNPs
env131 <- testPdp[[1]] %>% filter(response.id=='env_131')
summary(env131)
 ggplot(env131, aes(Longitude, yhat, color=response.id))+
  geom_line() +
  theme_bw()


#adding other response varables to see if this improves predictions. wecould say that it only has been tested 
#on 3 algorithms (GAMs, Xgboost which is an improved GBM and liner models) and user beware otherwise.

covariance_mod <- 
  boost_tree() %>%
  set_engine("xgboost") %>%
  set_mode("regression")

#or 
covariance_mod <- 
  rand_forest(trees = 100, mode = "regression") %>%
  set_engine("ranger", importance = "permutation") %>%
  set_mode("regression")

#rhats <- StackPredictions(yhats, data.test, covariance_mod = 'gam')

#calculate covariance predictions. Nick - I just split this function up as I imagine it would be handy.
#at this bit,
rhats_2 <- response_covariance(yhats, covariance_mod)

#stack everything together
finalPred <- stacked_preds(rhats_2, yhats)

#add similarity in predictions to look for SNP patterns.

#at some point it would be good to add interactions/Shapely here too. 

#---------------------------------------------------------------------------------------------------------------
