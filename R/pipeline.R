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
library(vegan)
library(ggrepel)
library(LEA) 

# load all function codes. This will disappear when we formally make this a function
source("./R/filterRareCommon.R")
source("./R/MrTidyModels.R")
source("./R/StackPredictions.R")
source("./R/devianceResids.R")
source("./R/filterRareCommon.R")
source("./R/MrtTidyPerf.R")
source("./R/mrvip.R")
source("./R/plot_vi.R")
source("./R/mrPdP.R")

#new interaction code
source("./R/mrInteractions.R")
source("./R/mrPlotInteractions.R") #not finding this for some reason - no idea why
source("./R/vintTidy.R")

#Nick C - this and the function below are the functions with tidy model code I made from your original
source("./R/stacked_preds.R") #this does the stacking
source("./R/response_covariance") #this should create the covatiance matrix

source("./R/readSnpsPed.R") #function for reading SNP data from plink .ped file
source("./R/ResistanceComponents.R") #function for generating resistance component data from resistance matrices

#----------------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------
#Example: Bobcat SNP data- from Plink file
#---------------------------------------------------------------------------------

snps <- readSnpsPed("bobcat.plink.ped", "bobcat.plink.map") #NAs in data and interpolated as the mode. 
X <- filterRareCommon (snps, lower=0.4, higher=0.75) #these are harsh

X <- X[-c(279:281),] #these didnt match

#Create resistance components for features using PCoA

#load resistance matrix generated using circuitscape. Add any matrix to a sub file in the working directory
#and provide that name to the function below. This will create a resistance component dataframe. Could add Mems here too

Y <- resist_components(data_resist,siteData, p_val=0.001,  filename = 'Bobcat_cs_matrices')

df1 <- mutate_all(Y, function(x) as.numeric(as.character(x))) #this is now in the function


#check for correlations and remove
df2 <- cor(df1) #find correlations
hc <-  findCorrelation(df2, cutoff=0.9) # put any value as a "cutoff". This is quite high
hc <-  sort(hc)

Y <-  df1[,-c(hc)] #reduced set of Y features. Labelled Y for simplicity

#---------------------------------------------------------------------------------
#Example: Bobcat FIV data  - read  viral data already curated
#---------------------------------------------------------------------------------

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

#for more efficent testing for interactions (more variables more interacting pairs)
Y <- FeaturesnoNA[c(1:3)]
#Optional: Filter rare/common SNPs or species. Retaining minor allelle frequncies >0.1 and removing common allelles (occur>0.9)
fData <- filterRareCommon (Responsedata, lower=0.4, higher=0.75) 
#for the snps data
fData <- filterRareCommon (snps, lower=0.4, higher=0.75) 

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
X <- read.csv('grid101.csv', row.names = NULL, head=T)
X[1:3] <- NULL
X[2:6] <- NULL; X[1:2] <- NULL #removing these for the moment.

X[X==2] <- 1 #'cough' this is terrible  - it would be easier to have each loci as two columns to make this properly 
 #a classification model (rather than dealing with nominal data but i need to rethink probability)

#remove NAs
X[X==NA] <- 0 #not optimal

#heavy filtering
Xsim <- filterRareCommon (X, lower=0.01, higher=0.999)  #this isnt working

Y <- read.csv('simple_sims_env.csv')

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
  rand_forest(trees = 100, mode = "classification") %>%
  set_engine("ranger", importance = c("impurity","impurity_corrected")) %>%
  set_mode("classification")

#model 1 with tuning. The user has to supply this too
model1tune<- rand_forest(
  mtry = tune(),
  trees = 100,
  min_n = tune()
) %>%
  set_mode("classification") %>%
  set_engine("ranger")

#Boosted model (xgboost). Not tuned either currently. This doesn't work with the small bobcat dataset. Much too small

model1 <- 
boost_tree() %>%
  set_engine("xgboost") %>%
  set_mode("classification")
  
#---------------------------------------------------------------------
#Perform the analysis
#---------------------------------------------------------------------  
#parallell processing
all_cores <- parallel::detectCores(logical = FALSE)

library(doParallel)
cl <- makePSOCKcluster(all_cores)
registerDoParallel(cl)
  
#models just using features/predictor variables.
yhats <- mrIMLpredicts(X=X,Y=Y, model1=model1, balance_data='no') ## in MrTidymodels. Balanced data= up updamples and down downsampled to create a balanced set

#-------------------------------------------------------------------
# Visualization for model tunning and performance
#-------------------------------------------------------------------
## GM add here

#we can now assess model performance from the best tuned model
ModelPerf <- mrIMLperformance(yhats, model1, X=X) # MCC is useful (higher numbers = better fit)

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


#-------------------------------------------------------------------
# Visualization individual and global feature importance
#-------------------------------------------------------------------
## GM add here

#we can look at variable importance. Co-infection data only has one feature so not much use there.
VI <- mrVip(yhats, Y=Y) #seem not to used all features for log regression for some reason
#plot model similarity

#plot variable importance

#for interpretation group features. Annoying but I cant see an easier way. Features are in alphabetical order
groupCov <- c(rep ("Host_characteristics", 1),rep("Urbanisation", 3), rep("Vegetation", 2), rep("Urbanisation",1), rep("Spatial", 2), 
              rep('Host_relatedness', 6),rep ("Host_characteristics", 1),rep("Vegetation", 2), rep("Urbanisation",1))  

#not that GLMs in particular wont produce coefficents for features that are strongly colinear and will drop them from the model.
#in this case group cov will have to be changed to reflect features included in the model. 

plot_vi(VI=VI,  X=X,Y=Y, modelPerf=ModelPerf, groupCov=groupCov, cutoff= 0.5)#note there are two plots here. 

#if you dont want/need to group covariates:

plot_vi(VI=VI,  X=X,Y=Y, modelPerf=ModelPerf, cutoff= 0.5) #mcc cutoff not working right

#First plot is overall importance and the second is individual SNP models.
#warning are about x axis labels so can ignore 

#plot partial dependencies. Strange results
testPdp <- mrPdP(yhats, X=X,Y=Y, Feature='Grassland') 

## make one plot for each 
## pre plot for each

#indiv SNPs
testPdp %>% 
  filter(response.id=="env_212")%>%
  ggplot(aes(Grassland , yhat, group = yhat.id))+
  geom_line()+
  theme_bw()

#calculate interactions  -this is qute slow and memory intensive

interactions <-mrInteractions(yhats, X, Y) 

mrPlot_interactions(Interact, X,Y, top_ranking = 3, top_response=3)

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
