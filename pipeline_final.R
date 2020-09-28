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
#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")

#BiocManager::install("LEA")
library(ape)
library(flashlight)
library(devtools)
install_github("mayer79/flashlight") 
#library(pbapply)

# load all function codes. This will disappear when we formally make this a function
source("./R/filterRareCommon.R")
source("./R/MrIMLpredicts.R")
source("./R/StackPredictions.R")
source("./R/devianceResids.R")
source("./R/filterRareCommon.R")
source("./R/MrIMLperformance.R")
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

source(("./R/mrFlashlight.R"))
source("./R/mrProfileplots.R")
#----------------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------
#Example: Bobcat SNP data- from Plink file
#---------------------------------------------------------------------------------

snps <- readSnpsPed("bobcat.plink.ped", "bobcat.plink.map.map") #NAs in data and interpolated as the mode. 
X <- filterRareCommon (snps, lower=0.4, higher=0.7) #these are harsh

X <- X[-c(279:281),] #these didnt match

#whatever individual X128 is it needs to go as it is a crazy outlier!

#filter correlated SNPS

df2 <- cor(X) #find correlations
hc <-  findCorrelation(df2, cutoff=0.5) # put any value as a "cutoff". This is quite high
hc <-  sort(hc)

# Labelled Y for simplicity

#Create resistance components for features using PCoA

#load resistance matrix generated using circuitscape. Add any matrix to a sub file in the working directory
#and provide that name to the function below. This will create a resistance component dataframe. Could add Mems here too

Y <- resist_components(data_resist, p_val=0.01,  filename = 'Bobcat_cs_matrices')

#check for correlations and remove
df2 <- cor(Y) #find correlations
hc <-  findCorrelation(df2, cutoff=0.9) # put any value as a "cutoff". This is quite high
hc <-  sort(hc)

Y <-  Y[,-c(hc)] #reduced set of Y features. Labelled Y for simplicity

#---------------------------------------------------------------------------------
#Example: Puma SNP data
#--------------------------------------------------------------------------------

snps <- readSnpsPed("ws.ped", "ws.map") #NAs in data and interpolated as the mode. 
X <- filterRareCommon (snps, lower=0.4, higher=0.75) #these are harsh

#---------------------------------------------------------------------------------
#Example: Bobcat FIV data  - read  viral data already curated
#---------------------------------------------------------------------------------

#load data for training
Responsedata <-read.csv("GF SNPs posBinary.csv")
rownames(Responsedata) <-Responsedata$X
Responsedata$X <- NULL
Responsedata[Responsedata==4] <-0 #fix random 4s in response data. 

#feature set. Note that samples must be rows.
Features <-read.csv("Landscape and host data.csv", row.names = 1, head=T)
#summary(Features)

# instead we can devlop a cut % for interpolation in the func_load_data. Seems to be a random row somehow added to the feature data.
#Features<-na.interpolation(Features, option = "spline")

# # or remove NAs from the feature/predictor data.
FeaturesnoNA<-Features[complete.cases(Features), ];str(Features) #dropping NAs

Y <- FeaturesnoNA #for simplicity

#for more efficent testing for interactions (more variables more interacting pairs)
Y <- FeaturesnoNA[c(1:3)]
#Optional: Filter rare/common SNPs or species. Retaining minor allelle frequncies >0.1 and removing common allelles (occur>0.9)
fData <- filterRareCommon (Responsedata, lower=0.4, higher=0.7) 
X <- fData #for simplicity when comparing

#that occur in less than 35% of individuals and > 75% of individuals
#fData <- rownames_to_column(fData, "Individual")#get individual id back
#user can specify any tidy model here. 

#---------------------------------------------------------------------------------
#coinfection test data # from MRFcov. This data can also be used in this pipeline replacing X/Y
#---------------------------------------------------------------------------------

X <- select(Bird.parasites, -scale.prop.zos) #response variables eg. SNPs, pathogens, species....
Y <- select(Bird.parasites, scale.prop.zos) # feature set

#---------------------------------------------------------------------------------
#Simulated data
#---------------------------------------------------------------------------------
X <- read.csv('grid101_binary.csv', row.names = NULL, head=T)
X[1:9] <- NULL
#heavy filtering
Xsim <- filterRareCommon (X, lower=0, higher=1)  #this isnt working for mutinomial

Y <- read.csv('simple_sims_env.csv')

#note that for multinomial models we do not have a great solution yet to calculate model deivance
#this means that model stacking currently is unavialable for this class of models.

#multinominal models arent working either...


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
#--------------------------------------------------------------------------------------------------------------------------------------
#creating the models

#-------------------------------------------------------------------
#Pre-training data visualization- checks for long tall distributions
#-------------------------------------------------------------------
## GM add here



#-------------------------------------------------------------------
#Set up the model
#-------------------------------------------------------------------

#user provides the models they would like for each component (model1)
 
model1 <- #model used to generate yhat
  # specify that the model is a random forest
  logistic_reg() %>%
  # select the engine/package that underlies the model
  set_engine("glm") %>%
  # choose either the continuous regression or binary classification mode
  set_mode("classification") #just for your response

model1 <- 
  rand_forest(trees = 100, mode = "classification") %>% #this should cope with multinomial data alreadf
  set_engine("ranger", importance = c("impurity","impurity_corrected")) %>%
  set_mode("classification")

#Boosted model (xgboost). Not tuned either currently. This doesn't work with the small bobcat dataset. Much too small

model1 <- 
boost_tree() %>%
  set_engine("xgboost") %>%
  set_mode("classification")

#for SVM need different tuning paramters. Currently tuning works best for tree-based algorithms.
  
#---------------------------------------------------------------------
#Perform the analysis
#---------------------------------------------------------------------  
#parallell processing

#models just using features/predictor variables.
yhats <- mrIMLpredicts(X=X,Y=Y, model1=model1, balance_data='no', model='classification', parallel = TRUE) ## in MrTidymodels. Balanced data= up updamples and down downsampled to create a balanced set. For regression there no has to be selected.

#-------------------------------------------------------------------
# Visualization for model tunning and performance
#-------------------------------------------------------------------
## GM add here

#we can now assess model performance from the best tuned model
 # MCC is useful (higher numbers = better fit)

yhats <- mrIMLpredicts(X=X,Y=Y, model1=model1, balance_data='no', model='classification', parallel = TRUE) ## in MrTidymodels. Balanced data= up updamples and down downsampled to create a balanced set. For regression there no has to be selected.

ModelPerf <- mrIMLperformance(yhats, model1, X=X, model='regression')
ModelPerf[[1]] #predictive performance for individual responses 
ModelPerf[[2]]#overall predictive performance. r2 for regression and MCC for classification


#-------------------------------------------------------------------
# Visualization individual and global feature importance
#-------------------------------------------------------------------
## GM add here

#we can look at variable importance. Co-infection data only has one feature so not much use there.
VI <- mrVip(yhats, Y=Y) 
#plot model similarity

#plot variable importance

#for interpretation group features. Annoying but I cant see an easier way. Features are in alphabetical order
groupCov <- c(rep ("Host_characteristics", 1),rep("Urbanisation", 3), rep("Vegetation", 2), rep("Urbanisation",1), rep("Spatial", 2), 
              rep('Host_relatedness', 6),rep ("Host_characteristics", 1),rep("Vegetation", 2), rep("Urbanisation",1))  

#not that GLMs in particular wont produce coefficents for features that are strongly colinear and will drop them from the model.
#in this case group cov will have to be changed to reflect features included in the model. 

plot_vi(VI=VI,  X=X,Y=Y, modelPerf=ModelPerf, groupCov=groupCov, cutoff= 0.5, plot.pca='no')#note there are two plots here. PCA is hard to read with > 50 response varianbles

#if you dont want/need to group covariates:

plot_vi(VI=VI,  X=X,Y=Y, modelPerf=ModelPerf, cutoff= 0, plot.pca='yes', model='regression') #mcc cutoff not working right

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

flashlightObj <- mrFlashlight(yhats, X, Y, response = "multi", model='regression')

#plots


plot(light_profile(flashlightObj, v = "bio_1", type = "ale"))
plot(light_profile(mfl, v = "bio_1", type = "ale"))
#this will plot all responses. See mrALEplots below

#another way to assess performance for each response. Lots of response make it hard to read
plot(light_performance(multifl), rotate_x = TRUE, fill = "orange") +
  
  labs(x = element_blank())


#plot prediction scatter for all responses. Gets busy with 
plot(light_scatter(flashlightObj, v = "Grassland", type = "predicted"))

#plots everything on one plot (partial dependency, ALE, scatter)
plot(light_effects(flashlightObj, v = "Grassland"), use = "all")


profileData_pd <- light_profile(flashlightObj, v = "bio_1") #partial dependencies
profileData_ale <- light_profile(flashlightObj, v = "bio_1", type = "ale") #acumulated local effects

#Plot global ALE plot. This the plots the smoothed average ALE value. Have to remove responses that don't respond.

mrProfileplot(profileData_pd , sdthresh =0.01)
mrProfileplot(profileData_ale , sdthresh =0.01)


#-------------------------------------------------------------------------------------------------

#calculate interactions  -this is qute slow and memory intensive

interactions <-mrInteractions(yhats, X, Y,  mod='regression') #this is computationally intensive so multicores are needed. If stopped prematurely - have to reload things

mrPlot_interactions(interactions, X,Y, top_ranking = 5, top_response=5)

save(interactions, 'Fitzpatrick2016interactions')

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
