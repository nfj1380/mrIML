
pacman::p_load('MRFcov', 'mrIML', 'tidyverse', 'future.apply','tidymodels', 'finetune',
               'themis', 'vip', 'flashlight', 'iml', 'vivid', 'igraph', 'ggnetwork', 'network',
                'gridExtra', 'xgboost', 'brulee', 'fastshap', 'tabnet', 'bonsai',
               'parsnip', 'cowplot', 'progress', 'hstats', 'geosphere')

source("./R/MrIMLpredicts.R")
source("./R/mrResponseStacks.R")
source("./R/mrCoOccur.R")
source("./R/mrCoOccurImp.R")
source("./R/mrCoOccurNet.R")
source("./R/mrBootstrap.R")
source("./R/mrPD_bootstrap.R")
source("./R/mrVI_bootstrap.R")
source("./R/mrCoOccurNet_bootstrap.R")
source("./R/mrCovar.R") 
source("./R/mrInteractionsSept23.R") 

#extra git package to download
#devtools::install_github("mayer79/hstats")

#---------------------------------------------------------------------------------
#coinfection test data # from MRFcov. This data can also be used in this pipeline replacing X/Y
#---------------------------------------------------------------------------------

Y <- dplyr::select(Bird.parasites, -scale.prop.zos) #response variables eg. SNPs, pathogens, species....
X <- dplyr::select(Bird.parasites, scale.prop.zos) # feature set

X1 <- Y %>%
  dplyr::select(sort(names(.)))

#need to make sure responses are arranged alphabetically

Y <- Y %>%
  dplyr::select(sort(names(.)))

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
Y <- filterRareCommon(Responsedata, lower=0.2, higher=0.8) %>% 
  dplyr::select(sort(names(.))) #0.1 for final analysis
#set up interactions
X1 <- Y 

############################
#Spatial data
############################

raw_spatial <- readRDS('coord_updated') #read data in
 
distance_matrix <- geosphere::distm(raw_spatial) #create a distance matrix

#single distance (0km by default)
mems <- spatialRF::mem(distance.matrix = distance_matrix)

#rank them and take the top 5. Could use the complte set too. 
mem.rank <- spatialRF::rank_spatial_predictors(
  distance.matrix = distance_matrix,
  spatial.predictors.df = mems,
  ranking.method = "moran"
)
#other random effects. Site and plate location
rand_eff <- readRDS('randonEff')
#can add them to X

#this will control for site effects too (ticks sampled at each site have one location recorded)
Xsp <- mem.rank$spatial.predictors.df[1:5]
Xcombined_re <- cbind(X, Xsp, rand_eff)
Xcombined <- cbind(X, Xsp);str(Xcombined)
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
################################
#Configuring the model
################################
#Can try different combinations of things here. With and without X1 (taxa assocations)
#With Xsp (spatial eigenvectors) and X (enviro/host variables) combined or sepparated
#compare performance across all data and a few differing algorithms

yhats_rf_noenviro <- mrIMLpredicts(X=NULL, Y=Y,
                          X1=X1,
                          spatial_data=raw_spatial,
                          Model=model_rf , #lm/xgb not working on avian malaria
                          balance_data='no',
                          mode='classification',
                          tune_grid_size=5,
                          seed = sample.int(1e8, 1),
                          morans=T,
                          prop=0.6, k=5, racing=T) #racing =F works better for small samples or low prev data

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


############################################################
#model performance and basic interpretation
yhats <- yhats_rf
ModelPerf <- mrIMLperformance(yhats, Model=model_rf, Y=Y, mode='classification')
ModelPerf[[2]]
m <- ModelPerf[[1]]

#check spatial Morans I - looking for taxa with siginificant MoranI (they are autocorrelated)

spat_list <- yhats_rf_noenviro  %>% purrr::map(pluck("moran_p")) 

#check the list and identify taxa with signal. The tick microbiome is
spat_data <- do.call(rbind,spat_list )

# Run bootraps -uses pdps as they are more flexible and easier to interpret. But have pdp issues if features are correlated

bs_tick <- mrBootstrap(yhats=yhats, Y=Y, num_bootstrap = 20, alpha = 0.05, ice=F) #make sure ice=F at this stage
bs_malaria <- mrBootstrap(yhats=yhats,Y=Y, num_bootstrap = 5, alpha = 0.05, ice=F) 


#plot bootstrap importance
bs_impVIa <- mrVI_bootstrap(mrBootstrap_obj=bs_malaria, ModelPerf=ModelPerf, 
                           threshold=0.8,  X=X, Y=Y, global_top_var=10,
                           local_top_var=5)
vi_obj <- bs_impVIa[[1]]#data for posterity
bs_impVIa[[2]]#combined plot

#create bootstrapped pdps. Global_top5 allows just to plot the # most important features

pds <- mrPD_bootstrap(mrBootstrap_obj=bs_malaria, vi_obj=bs_impVIa, X, Y,
                      target='Plas', global_top_var=5)
pd_list <- pds[[1]] #data
pds[[2]]#plot 

###########################################################################
#co-occurence and interactions
###########################################################################

#still updating this function
int_ <- mrInteractions(yhats, X, Y, num_bootstrap=100,
                           feature = 'Hkillangoi', top.int=10)

int_[[1]] # overall plot
int_[[2]] # individual plot for the response of choice 
int_[[3]] #two way plot

#can also look at three-way interactions for taxa of interest.
#Note these calues are not bootrapped 

#extract model
model_fit <- yhats[[2]]$mod1_k %>% extract_fit_parsnip() #second taxa

var_names <- names(yhats[[2]]$data)[-1] #second taxa



#pred function needs to be supplied
pred_fun <- function(m, dat) {
  predict(
    m, dat[, colnames(yhats[[2]]$data)[-1], drop = FALSE],
    type = "prob"
  )$`.pred_1`
}

s <- hstats(model_fit, v = names(yhats[[2]]$data_train)[-1], #make sure taxa match (e.g. ,1)
            X = yhats[[2]]$data_train, pred_fun = pred_fun, n_max = 300, 
            pairwise_m = 10)

summary(s) #provides summary output

plot(s, which = 1:3, normalize = F,
     squared = F, facet_scales = "free_y", ncol = 1)

#plot some interactions.Still a work in progress

####NB doesnt work
#pd <- partial_dep(model_fit, v = c("scale.prop.zos", "Plas"), X = yhats[[2]]$data_train)

fl <- mrFlashlight(yhats, X=cbind(X, Y), Y=Y, response = "single", index=2, mode='classification') #crate a wrapper than just constructs lots of dingles

plot(light_profile2d(fl, c("scale.prop.zos", "Plas")))+theme_bw()
  #add a rug plot
#convert to a iml object

iml_convert_list <- MrIMLconverts(yhats, X=cbind(X,Y),  mode='classification')

#still in progress. Old function needs to be put into iml package
test <- mrProfile2D(iml_convert_list, featureA='scale.prop.zos', featureB='Plas',
                        mode='classification', grid.size=30, method = "ale")


plot(light_profile2d(fl, c("Plas", "Hzosteropis")))
#create association network

assoc_net<- mrCoOccurNet_bootstrap (mrPD_obj=pds , Y=Y)

####################
#Example Plot
####################

#we propose the following rule of thumb for associations. 0.05 threshold (mean strength).
#this could be tailored though  -anyone got any better ideas? Ideally we'd have bootrap values
#on each edge, but it isnt possible with the current approach. Could us minium too. 0.01 works ok

assoc_net_filtered <-  assoc_net %>% 
  filter(mean_strength > 0.05)


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

fl <- mrFlashlight(yhats, X=cbind(X, Y), Y=Y, response = "single", index=2, mode='classification') #crate a wrapper than just constructs lots of dingles

treeA <- light_global_surrogate(fl)
treeA$data
plot(treeA)

#many other ways to interpret individual responses

###########################################################################
#Exploring covariates
###########################################################################

#explore the environmental/host covariates in more detail.
#factor box plots need to be fixed
covar <- mr_Covar(yhats, X=X, Y=Y, var='scale.prop.zos', sdthresh =0.01) #sdthrsh just plots taxa responding the most.
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

yhats_rf_stacks <- mrResponseStacks(yhats_rf_noX1, taxa=Y, Model=model_rf_reg, seed=123)

ModelPerf <- mrIMLperformance(yhats=yhats_rf_stacks , Model=model_rf, Y=Y, mode='regression')

VI <- mrCoOccur(yhats=yhats, taxa=Y, X=X) 

#calculate deviance
dList <- yhats %>% purrr::map(pluck("deviance"))
deviance_df <- data.frame(do.call(cbind, dList))
deviance_df_noinf<- deviance_df %>% mutate_all(~ ifelse(. == -Inf, -3, .))
deviance_df_noinf<- deviance_df_noinf %>% mutate_all(~ ifelse(. == Inf, 3, .))
# glimpse(deviance_df_noinf)

final_deviance <- deviance_df_noinf

fl <- mrFlashlight(yhats_rf_stacks , X=X1, Y=final_deviance, response = "single", index=1, mode='regression')
a = (light_scatter(fl, v = "Hkillangoi", type = "predicted"))
a +geom_boxplot()
plot(light_global_surrogate(fl))

plot(light_profile2d(fl, c("Plas", "Hzosteropis")))

#plot(light_profile(fl, v = "Hkillangoi", type = "ale"))
#below doesnt work
# profileData_ale <- light_profile(fl, v = "Hkillangoi", type = "ale")
# test <- profileData_ale$data
# mrProfileplot(profileData_ale  , sdthresh =0) 
