
pacman::p_load('MRFcov', 'mrIML', 'tidyverse', 'future.apply','tidymodels', 'finetune',
               'themis', 'vip', 'flashlight', 'iml', 'vivid', 'igraph', 'ggnetwork', 'network',
               'gridExtra', 'xgboost', 'brulee', 'fastshap', 'tabnet', 'bonsai',
               'parsnip', 'cowplot', 'progress', 'hstats')

source("./R/MrIMLpredicts.R")
source("~/MrIML/mrIML/R/mrResponseStacks.R")
source("./R/mrCoOccur.R")

source("./R/mrCoOccurNet.R")
source("./R/mrBootstrap.R")
source("./R/mrPD_bootstrap.R")
source("./R/mrVI_bootstrap.R")
source("./R/mrCoOccurNet_bootstrap.R")
source("./R/mrCovar.R") 

source("./R/mrResponseStacks.R")
#extra git package to download
#devtools::install_github("mayer79/hstats")

#---------------------------------------------------------------------------------
#coinfection test data # from MRFcov. This data can also be used in this pipeline replacing X/Y
#---------------------------------------------------------------------------------

Y <- select(Bird.parasites, -scale.prop.zos) #response variables eg. SNPs, pathogens, species....
X <- select(Bird.parasites, scale.prop.zos) # feature set

X1 <- Y %>%
  select(sort(names(.)))

#need to make sure responses are arranged alphabetically

Y <- Y %>%
  select(sort(names(.)))

# Tick microbiome data
#---------------------------------------------------------------------------------

FeaturesnoNA <- read.csv('tick_covar.csv') %>% glimpse()
#data cleaning
FeaturesnoNA$X <- NULL

#for more efficent testing for interactions (more variables more interacting pairs)
X <- FeaturesnoNA
#DataExplorer::introduce(X)

#responses

Responsedata<- read.csv('Tick_microbe_asvs.csv') %>% glimpse()

#dont include rare and very common species
Y <- filterRareCommon(Responsedata, lower=0.3, higher=0.7) %>% 
  select(sort(names(.))) #0.1 for final analysis
#set up interactions
X1 <- Y 

#---------------------------------------------------------------------------------

#model setup
model_rf <- 
  rand_forest(trees = 100, mode = "classification", mtry = tune(), min_n = tune()) %>% #100 trees are set for brevity. Aim to start with 1000
  set_engine("randomForest")

model_lm <- #model used to generate yhat
  logistic_reg() %>%
  set_engine("glm") %>%
  set_mode("classification") #just for your response

#doesn't work  - maybe catboost? or H2o
model_xgb <- 
  boost_tree(
    mtry = tune(), trees = 100, tree_depth = tune(), 
    learn_rate = tune(), min_n = tune(), loss_reduction =tune()) %>% 
  set_mode("classification")

model_dnn <-
  mlp(
    hidden_units = 10, # 10 to start with
    dropout      = tune(),
    epochs       = tune(),
    learn_rate   = tune(),
    activation   = "elu"
  ) %>% 
  set_engine("brulee") |>
  set_mode("classification")


cl <- parallel::makeCluster(5)
plan(cluster, workers=cl)
#
#not working for xgb or lms  - dummy problem
yhats_dnn <- mrIMLpredicts(X=X,
                           Y=Y,#X1=X1,
                           Model=model_dnn  , #lm doesnt work on co-occurence models
                           balance_data='no',
                           mode='classification',
                           tune_grid_size=5,
                           seed = sample.int(1e8, 1),
                           prop=0.7, k=5, racing=T)

yhats_xgb <- mrIMLpredicts(X=X,
                           Y=Y, X1=X1,
                           Model=model_xgb , #lm/xgb not working on avian malaria
                           balance_data='down',
                           mode='classification',
                           tune_grid_size=5,
                           seed = sample.int(1e8, 1),
                           prop=0.7, racing=T)

yhats_rf <- mrIMLpredicts(X=X,
                          Y=Y, X1=X1,
                          Model=model_rf , #lm/xgb not working on avian malaria
                          balance_data='no',
                          mode='classification',
                          tune_grid_size=5,
                          seed = sample.int(1e8, 1),
                          prop=0.6, k=5, racing=T) #racing =F works better for small samples or low prev data

############################################################
#model performance and basic interpretation

ModelPerf <- mrIMLperformance(yhats_rf, Model=model_rf, Y=Y, mode='classification')
m <- ModelPerf[[1]]

#check importance of complete model. 

# Run bootraps -uses pdps as they are more flexible and easier to interpret. But have pdp issues if features are correlated

bs_tick <- mrBootstrap(yhats=yhats_rf, Y=Y, num_bootstrap = 20, alpha = 0.05, ice=F) #make sure ice=F at this stage
bs_malaria <- mrBootstrap(yhats=yhats,Y=Y, num_bootstrap = 5, alpha = 0.05, ice=F) 


#plot bootstrap importance
bs_impVIa <- mrVI_bootstrap(mrBootstrap_obj=bs_tick, ModelPerf=ModelPerf, 
                            threshold=0.9,  X=X, Y=Y, global_top_var=10,
                            local_top_var=5)
vi_obj <- bs_impVIa[[1]]#data for posterity
bs_impVIa[[2]]#combined plot

#create bootstrapped pdps. Global_top5 allows just to plot the # most important features

pds <- mrPD_bootstrap(mrBootstrap_obj=bs_tick, vi_obj=bs_impVIa, X, Y,
                      target='Bacteroides', global_top_var=5)
pds[[1]] #data
pds[[2]]#plot

###########################################################################
#co-occurence and interactions
###########################################################################

#still updating this function
int_test <- mrInteractions(yhats,  ModelPerf=ModelPerf, X, Y, num_bootstrap=10)

#still a work in progress wit boots but works 
assoc_net<- mrCoOccurNet_bootstrap (mrPD_obj=pds , Y=Y)

####################
#Example Plot
####################

#we propose the following rule of thumb for associations. 0.05 threshold (mean strength).
#this could be tailored though  -anyone got any better ideas? Ideally we'd have bootrap values
#on each edge, but it isnt possible with the current approach. Could us minium too. 0.01 works ok

assoc_net_filtered <-  assoc_net %>% 
  filter(mean_strength > 0.01)


#convert to igraph
g <- graph_from_data_frame(assoc_net_filtered, directed=TRUE, vertices=names(Y)) #matching Y data

#need to fix this - reformat the dataframe
#plot(g)

E(g)$Value <- assoc_net_filtered$mean_strength###chnge this as needed
E(g)$Color <- ifelse(assoc_net_filtered$direction == "negative", "blue", "red")

# Obtain NMDS layout coordinates
#nmds_layout <- layout_with_mds(g, dim = 2) #this isnt great on many

# Convert the igraph object to a ggplot object with NMDS layout
gg <- ggnetwork(g)

# Plot the graph
ggplot(gg, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(aes(color = Color, linewidth = (Value)), curvature = 0.2,
             arrow = arrow(length = unit(5, "pt"), type = "closed")) + #makes arrows bigger
  geom_nodes(color = "gray", size = degree(g, mode = "out")/2)+#, fill='black', stroke =2) +
  scale_color_identity() +
  theme_void() +
  theme(legend.position = "none")  +
  geom_nodelabel_repel(aes(label = name),
                       box.padding = unit(0.5, "lines"),
                       data = gg,
                       size=2,
                       segment.colour = "black",
                       colour = "white", fill = "grey36")

###########################################################################
#individual taxa level
###########################################################################

#str(Y)

fl <- mrFlashlight(yhats_rf, X=cbind(X, Y), Y=Y, response = "single", index=5, mode='classification') #crate a wrapper than just constructs lots of dingles

treeA <- light_global_surrogate(fl)
treeA$data
plot(treeA)

#many other ways to interpret individual responses

#vivid. Doesnt work yet.
mrIMLconverts_list <- MrIMLconverts(yhats_rf, cbind(Y,X), mode='classification')

#got this on the other computer
vivi_test <- vivi(fit =  extract_fit_parsnip(yhats_rf1[[36]]$mod1_k),
                  data = yhats_rf1[[36]]$data,
                  response = "class",
                  reorder = FALSE,
                  normalized = FALSE)
###########################################################################
#Exploring covariates
###########################################################################

#explore the environmental/host covariates in more detail. No bootraps for these yet
#factor box plots need to be fixed
covar1 <- mr_Covar(yhats_rf, X=X, Y=Y, var='max_snow_depth', sdthresh =0.01) #sdthrsh just plots taxa responding the most.
#must add a warning her eif threshold is too high

#bottom plots shows the overall change in probabilities across each unit of X.

#################################
#stacked approach
#################################

#the direct joint approach above adding all taxa + environment to each model works well 
#for smaller datasets < 1000 taxa (or so still testing)

#must run the model this time without X1 ()
yhats_rf_noX1 <- mrIMLpredicts(X=X,
                               Y=Y,
                               Model=model_rf , #lm/xgb not working on avian malaria
                               balance_data='no',
                               mode='classification',
                               tune_grid_size=5,
                               seed = sample.int(1e8, 1),
                               prop=0.6, k=5, racing=T) #


#must use a regression model to model the deviance from the first 

model_rf_reg <- 
  rand_forest(trees = 100, mode = "regression", mtry = tune(), min_n = tune()) %>% #100 trees are set for brevity. Aim to start with 1000
  set_engine("randomForest")

yhats_rf_stacks <- mrResponseStacks(yhats=yhats_rf_noX1, taxa=Y, Model=model_rf_reg, seed=123)

ModelPerf <- mrIMLperformance(yhats_rf_stacks , Model=model_rf, Y=Y, mode='regression')

VI <- mrStacksVI(yhats_stacks=yhats_rf_stacks,  Y=Y,  ice = FALSE, local_top_var=5) 
VI[[2]]

fl <- mrFlashlight(yhats_rf_stacks , X=taxa, Y=final_deviance, response = "single", index=1, mode='regression')
a = (light_scatter(fl, v = "Hkillangoi", type = "predicted"))
a +geom_boxplot()
plot(light_global_surrogate(fl))

plot(light_profile2d(fl, c("Plas", "Hzosteropis")))

#plot(light_profile(fl, v = "Hkillangoi", type = "ale"))
#below doesnt work
# profileData_ale <- light_profile(fl, v = "Hkillangoi", type = "ale")
# test <- profileData_ale$data
# mrProfileplot(profileData_ale  , sdthresh =0) 
