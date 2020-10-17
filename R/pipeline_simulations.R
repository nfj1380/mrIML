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
library(ape)
library(flashlight)
library(devtools)
install_github("mayer79/flashlight") 
#library(pbapply)

# load all function codes. This will disappear when we formally make this a function
source("./R/filterRareCommon.R")
source("./R/mrIMLpredicts.R")
source("./R/StackPredictions.R")
source("./R/devianceResids.R")
source("./R/filterRareCommon.R")
source("./R/mrIMLperformance.R")
source("./R/mrvip.R")
source("./R/plot_vi.R")

#new interaction code
source("./R/mrInteractions.R")
source("./R/mrPlotInteractions.R") #not finding this for some reason - no idea why
source("./R/vintTidy.R")

#Nick C - this and the function below are the functions with tidy model code I made from your original
source("./R/stacked_preds.R") #this does the stacking
source("./R/response_covariance") #this should create the covatiance matrix #not working?

source("./R/readSnpsPed.R") #function for reading SNP data from plink .ped file
source("./R/ResistanceComponents.R") #function for generating resistance component data from resistance matrices

source(("./R/MrFlashLight.R"))
source("./R/MrProfilePlots.R")
#----------------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------
#Simulated data
#---------------------------------------------------------------------------------
X <- read.csv('grid101_binary.csv', row.names = NULL, head=T)
X[1:9] <- NULL

#Replace NAs with major allele
Mode <- function(x, na.rm = FALSE) {
  if(na.rm){
    x = x[!is.na(x)]
  }
  
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}

for(i in 1:ncol(X)){
  X[i] <- na.replace(X[i], Mode(X[i], na.rm=TRUE)) #replace NAs with the major allele at a given locus
}




Y <- resist_components( p_val=0.01,  filename = 'Sim_matrices')


#note that for multinomial models we do not have a great solution yet to calculate model deivance
#this means that model stacking currently is unavialable for this class of models.

#multinominal models arent working either...


###############################################################################
#--------------------------------------------------------------------------------------------------------------------------------------
#creating the models

#-------------------------------------------------------------------
#Pre-training data visualization- checks for long tall distributions
#-------------------------------------------------------------------
## GM add here



#-------------------------------------------------------------------
#Set up the model
#-------------------------------------------------------------------

#Try random forest first, then XGR boost and log regression

model1 <- 
  rand_forest(trees = 100, mode = "classification") %>% #this should cope with multinomial data alreadf
  set_engine("ranger", importance = c("impurity","impurity_corrected")) %>%
  set_mode("classification")

model2 <- 
  boost_tree() %>%
  set_engine("xgboost") %>%
  set_mode("classification")

model3 <- #model used to generate yhat
  logistic_reg() %>%
  set_engine("glm") %>%
  set_mode("classification") #just for your response


#---------------------------------------------------------------------
#Perform the analysis
#---------------------------------------------------------------------  

#models just using features/predictor variables.
yhats_rf <- mrIMLpredicts(X=X,Y=Y, model1=model1, balance_data='no', model='classification', parallel = TRUE) ## in MrTidymodels. Balanced data= up updamples and down downsampled to create a balanced set. For regression there no has to be selected.
save(yhats_rf, file='./rf_model_sim_matrices')
ModelPerf <- mrIMLperformance(yhats_rf, model1, X=X, model='classification')
ModelPerf[[2]]
rm(yhats_rf)

yhats_xg <- mrIMLpredicts(X=X,Y=Y, model1=model2, balance_data='no', model='classification', parallel = TRUE) ## in MrTidymodels. Balanced data= up updamples and down downsampled to create a balanced set. For regression there no has to be selected.
save(yhats_xg, file='./xg_model_sim_matrices')
ModelPerf <- mrIMLperformance(yhats_xg, model2, X=X, model='classification')
ModelPerf[[2]]
rm(yhats_xg)

yhats_lr <- mrIMLpredicts(X=X,Y=Y, model1=model3, balance_data='no', model='classification', parallel = TRUE) ## in MrTidymodels. Balanced data= up updamples and down downsampled to create a balanced set. For regression there no has to be selected.
save(yhats_lr, file='./lr_model_sim_matrices')
ModelPerf <- mrIMLperformance(yhats_lr, model3, X=X, model='classification')
ModelPerf[[2]]
rm(yhats_lr)

#-------------------------------------------------------------------
# Visualization for model tunning and performance
#-------------------------------------------------------------------


load('rf_model_sim')


ModelPerf <- mrIMLperformance(yhats, model1, X=X, model='classification')
ModelPerf[[1]] #predictive performance for individual responses 
ModelPerf[[2]]#overall predictive performance. r2 for regression and MCC for classification

modelPerfD <- ModelPerf[[1]]

save(modelPerfD, file = 'simRF_performance')
#-------------------------------------------------------------------
# Visualization individual and global feature importance
#-------------------------------------------------------------------
## GM add here

#we can look at variable importance. Co-infection data only has one feature so not much use there.
VI <- mrVip(yhats, Y=Y) 
#plot model similarity

#plot variable importance

#not that GLMs in particular wont produce coefficents for features that are strongly colinear and will drop them from the model.
#in this case group cov will have to be changed to reflect features included in the model. 

plot_vi(VI=VI,  X=X,Y=Y, modelPerf=ModelPerf, cutoff= 0.6, plot.pca='no')#note there are two plots here. PCA is hard to read with > 50 response varianbles

#if you dont want/need to group covariates:

#First plot is overall importance, the second a pca showing responses with similar importance scores and the third is individual SNP models.
#warning are about x axis labels so can ignore 


# can build a flaslight object for individual responses 

fl <- mrFlashlight(yhats, X, Y, response = "single", index=1, model='classification')

plot(light_performance(fl), fill = "orange", rotate_x = TRUE) +
  
  labs(x = element_blank()) +
  
  theme(axis.text.x = element_text(size = 8))

str(fl)

#plot(light_breakdown(fl , new_obs = cbind(X, Y)[1, ]),by = X, v=Y) #prints all responses - need to fix but could be quite handy.

#int <- light_interaction(fl, pairwise=TRUE) #not working, but possible!

#------------------------------------------------------------
#Multiple response

flashlightObj <- mrFlashlight(yhats, X, Y, response = "multi", model='classification')

#plots

plot(light_profile(flashlightObj, v = "simple", type = "ale"))
plot(light_profile(mfl, v = "bio_1", type = "ale"))
#this will plot all responses. See mrALEplots below

#another way to assess performance for each response. Lots of response make it hard to read
plot(light_performance(multifl), rotate_x = TRUE, fill = "orange") +
  
  labs(x = element_blank())


#plot prediction scatter for all responses. Gets busy with 
plot(light_scatter(flashlightObj, v = "Grassland", type = "predicted"))

#plots everything on one plot (partial dependency, ALE, scatter)
plot(light_effects(flashlightObj, v = "Grassland"), use = "all")


profileData_pd <- light_profile(flashlightObj, v = "simple") #partial dependencies
profileData_ale <- light_profile(flashlightObj, v = "simple", type = "ale") #acumulated local effects

#Plot global ALE plot. This the plots the smoothed average ALE value. Have to remove responses that don't respond.

mrProfileplot(profileData_pd , sdthresh =0.07)
mrProfileplot(profileData_ale , sdthresh =0.07)


#-------------------------------------------------------------------------------------------------

#calculate interactions  -this is qute slow and memory intensive

interactions <-mrInteractions(yhats, X, Y,  mod='classification') #this is computationally intensive so multicores are needed. If stopped prematurely - have to reload things

mrPlot_interactions(interactions, X,Y, top_ranking = 5, top_response=5)

save(interactions, file='Fitzpatrick2016interactions')

#-------------------------------------------------------------------------------------------------


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

#Extra importance plots:

#another way to look at it
ModelPerf_p<-do.call(rbind.data.frame, ModelPerf)

## perfromance by outcome
ModelPerf_p%>%
  drop_na()%>%
  ggplot(aes(sensitivity, specificity, shape=response, colour=response  , fill=response)) +
  geom_smooth(method="lm") +
  geom_point(size=3)

## also
ModelPerf_p%>%
  drop_na()%>%
  ggplot(aes(mcc, specificity, shape=response, colour=response, fill=response   )) +
  geom_smooth(method="lm") +
  geom_point(size=3)

### Not sure we want to have this
ModelPerf_p[3] #summary mcc for all response variables
