#####################################################################
#PART 1: Microbiome simulation
#
# Contact:
# Kazuhiro Takemoto (takemoto@bio.kyutech.ac.jp)
#####################################################################

pacman::p_load('seqtime','SpiecEasi', 'igraph', 'mrIML', 'ppcor', 'PRROC')

pacman::p_load('MRFcov', 'mrIML', 'tidyverse', 'future.apply','tidymodels', 'finetune',
               'themis', 'vip', 'flashlight', 'iml', 'vivid', 'igraph', 'ggnetwork', 'network',
               'gridExtra', 'xgboost', 'brulee', 'fastshap', 'tabnet', 'bonsai',
               'parsnip', 'cowplot', 'progress', 'hstats', 'spatialRF','seqtime','SpiecEasi',
               'ppcor', 'PRROC')

# library(seqtime) # https://github.com/hallucigenia-sparsa/seqtime
# library(SpiecEasi) # https://github.com/zdk123/SpiecEasi
#install_github("zdk123/SpiecEasi")
source("generateM_specific_type.R")

# library(devtools)  
# install_github("hallucigenia-sparsa/seqtime")  
# library(seqtime) 

nn <- 40 # network size (n)
k_ave <- 4 # average degree (<k>). Could base this off real world nets? Try 2 and 4

## Generate an interaction matrix (Mij)
obj <- generateM_specific_type(nn,k_ave,type.network="random",type.interact="compt",interact.str.max=0.7,mix.compt.ratio=0.5)
#max interaction strength reduced from 0.8 to 0.5 for mutualism sims to ensure convergence
# @param nn number of nodes
# @param k_ave average degree (number of edges per node)
# @param type.network network structure
#               random: random networks
#                   sf: scale-free networks
#                   ws: small-world networks
#                   bipar: random bipartite networks
# @param type.interact interaction type
#               random: random
#               mutual: mutalism (+/+ interaction)
#                compt: competition (-/- interaction)
#                   pp: predator-prey (+/- interaction) - not suitable here
#                  mix: mixture of mutualism and competition. Mix 1 and 2 are the difficult ones
#                 mix2: mixture of competitive and antagonistic interactions
# @param interact.str.max maximum interaction strength
# @param mix.compt.ratio the ratio of competitive interactions to all intereactions (this parameter is only used for type.interact="mix" or ="mix2")


# adjacency matrix of network (Aij). Will compare this

network_real <- obj[[1]] #MRIML couldnt work out competition network. Mutualist doesnt work
# interaction matrix (Mij) for the GLV model
M <- obj[[2]]

## plot population dynamics
y <- rpois(nn,lambda=100)
r <- runif(nn)
res <- glv(nn, M, r, y)
tsplot(10*res[,20:1000],time.given =T)

#use defaults to geneerate dataset friom the glv model (300 samples and 150)
data <- as.data.frame(generateDataSet(300, M, count = nn*100, mode = 4)) #mode follows from Hirano and Takemoto 2019

#create a dataframe?
Mdf <-as.data.frame( obj[[2]])
colnames(Mdf) <- row.names(Mdf)

#convert to igraph
#diag(M) <- 0
diag(Mdf) <- 0
Mdf_mat <- as.matrix(Mdf) #to make it a matrix once again
g <- graph.adjacency(Mdf_mat , mode = "undirected", weighted = TRUE)#matching Y data

# Set edge colors based on positive (red) and negative (blue) weights
edge_colors <- ifelse(E(g)$weight > 0, "red", "blue")

# Set edge widths based on the absolute values of edge weights
edge_widths <- abs(E(g)$weight) * 5  # Adjust the scaling factor as needed

plot(
  g,
  edge.color = edge_colors,
  edge.width = edge_widths
)
gg <- ggnetwork(g)

gg2 <- gg %>% 
  mutate(edge_widths= ifelse(is.na(weight), 0, abs(weight))) %>% 
  cbind(edge_color = edge_colors)

# Plot the graph
ggplot(gg2, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(aes(color = edge_color, linewidth = (edge_widths)), curvature = 0.2,
             arrow = arrow(length = unit(5, "pt"), type = "closed")) + #makes arrows bigger
  geom_nodes(color = "gray", size = degree(g)/2)+
  scale_color_identity() +
  theme_void() +
  theme(legend.position = "none")  +
  geom_nodelabel_repel(aes(label = name),
                       box.padding = unit(0.5, "lines"),
                       data = gg,
                       size=2,
                       segment.colour = "black",
                       colour = "white", fill = "grey36")

#################################################
## Generate pa and relative abundance dataset
#################################################

# relative abundance
data_relative <- as.data.frame(t(data) / apply(data,2,sum))

 threshold <- 0.01 #threshold to say if a species is present < 1% of a sample =absent
 
# Transform the data to make presence absence
data_pa  <- data_relative %>%
  mutate(across(everything(), ~ ifelse(. <= threshold, 0, 1)))

##############################################################
#MrIML analysis
##############################################################
#mutualistic networks are the problem 
Y <- filterRareCommon(data_pa, lower=0.1, higher=0.9)  #9 taxa lost. Mostly had high abundance at each site

glimpse(Y)

#set up interactions
X1 <- Y 

model_rf <- 
  rand_forest(trees = 100, mode = "classification", mtry = tune(), min_n = tune()) %>% #100 trees are set for brevity. Aim to start with 1000
  set_engine("randomForest")

cl <- parallel::makeCluster(5)
plan(cluster, workers=cl)

yhats_rf_sim <- mrIMLpredicts(Y=Y, X=NULL,
                              X1=X1,
                              Model=model_rf , #lm/xgb not working on avian malaria
                              balance_data='no',
                              mode='classification',
                              tune_grid_size=5,
                              seed = sample.int(1e8, 1),
                              prop=0.7, k=5, racing=TRUE) #racing was having issues

ModelPerf <- mrIMLperformance(yhats_rf_sim , Model=model_rf, Y=Y, mode='classification')
m <- ModelPerf[[1]]
#all the species predicted poorly are not connected in the network for the competition

bs_sim <- mrBootstrap(yhats=yhats_rf_sim, Y=Y, num_bootstrap = 20, alpha = 0.05, ice=F)

bs_impVIa <- mrVI_bootstrap(mrBootstrap_obj=bs_sim, ModelPerf=ModelPerf, 
                            threshold=0.6,  X=X1, Y=Y, global_top_var=5,
                            local_top_var=5)
vi_obj <- bs_impVIa[[1]]#data for posterity
bs_impVIa[[2]]#combined plot

#this matches up with the real network

#create bootstrapped pdps. Global_top5 allows just to plot the # most important features

pds <- mrPD_bootstrap(mrBootstrap_obj=bs_sim, vi_obj=bs_impVIa, X=X1, Y=Y,
                      target='sp17', global_top_var=5)
pds[[1]] #data
pds[[2]]#plot

############################################
#Network
##############################################
assoc_net<- mrCoOccurNet_bootstrap (mrPD_obj=pds , Y=Y)

assoc_net_filtered <-  assoc_net %>% 
  filter(mean_strength > 0.05) #seems like the sweetspot? 

#to make everthing comparable
assoc_net_trans <- assoc_net_filtered  %>%
  mutate(mean_strength_dir = ifelse(direction == "negative", -mean_strength, mean_strength))


#convert to igraph
g <- graph_from_data_frame(assoc_net_trans , directed=FALSE, vertices=names(Y)) #matching Y data

edge_colors <- ifelse(E(g)$mean_strength_dir > 0, "red", "blue")

# Set edge widths based on the absolute values of edge weights
edge_widths <- abs(E(g)$mean_strength) * 5

plot(
  g,
  edge.color = edge_colors,
  edge.width = edge_widths
)

#need to fix this - reformat the dataframe
#plot(g)

E(g)$Value <- assoc_net_filtered$mean_strength###chnge this as needed
E(g)$Color <- ifelse(assoc_net_filtered$direction == "negative", "blue", "red")
E(g)$Value2 <- transformed_strength$mean_strength_dir #


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


###################################################
#How well does MrIML do?
###################################################

mrIML_matrix <- get.adjacency(g, attr = "mean_strength_dir", sparse = FALSE)
mrIML_matrix

#need to make sure everything matches for each network

# Remove 'sp' prefix from row names
new_row_names <- sub("^sp", "", rownames(mrIML_matrix))

# Set the new row names
rownames(mrIML_matrix) <- new_row_names

matching_indices <- match(rownames(mrIML_matrix), rownames(Mdf))
adj_M_reordered <- M[matching_indices, matching_indices]

#dist_matrix1 <- 1 - adj_M_reordered
dist_matrix1 <- 1 - adj_M_reordered #something is wronghere
dist_matrix_mrIML <- 1 - mrIML_matrix

#mantel test to see if these adjacency matrices are different or not
mantel_result <- vegan::mantel(dist_matrix1, dist_matrix_mrIML, method = "pearson")
mantel_result #not bad for a first attempt 0.52 and significant

###################################################
## Evaluating evaluating co-occurrence network performance
#using Mantel tests for other co-occurence model types
#Focus on SPIEC-EASI and SPacc
###################################################


# based on SPIEC-EASI
network_pred_spiec <- spiec.easi(t(data),method='mb')
# network_pred_spiec_mat <- as.matrix(getOptMerge(network_pred_spiec))
# dist_matrix_spiec <- 1 - network_pred_spiec_mat

network_pred_spiec_mat <- symBeta(getOptBeta(network_pred_spiec), mode='lower')

dissimilarity_matrix <- -1 / abs(network_pred_spiec_mat)
dissimilarity_matrix[dissimilarity_matrix == -Inf] <- 0

dist_matrix_spiec <- 1 - dissimilarity_matrix

#plot spiec
g <- graph.adjacency(network_pred_spiec_mat, mode = "undirected", weighted = TRUE)#matching Y data

# Set edge colors based on positive (red) and negative (blue) weights
edge_colors <- ifelse(E(g)$weight > 0, "red", "blue")

# Set edge widths based on the absolute values of edge weights
edge_widths <- abs(E(g)$weight) * 5  # Adjust the scaling factor as needed

plot(
  g,
  edge.color = edge_colors,
  edge.width = edge_widths
)
#on complete 'real' dataset

dist_matrix1a <- 1-M
diag(dist_matrix1a) <- 0 #maybe not needed? Doesnt change Mantel's
#SPIEC-EASI vs real
mantel_result1 <- vegan::mantel(dist_matrix1a, dist_matrix_spiec, method = "pearson")

#SPACC

network_pred_sparcc <- abs(sparcc(t(data))$Cor)

threshold_spacc <- 0.01 #keep this threshold for now
network_pred_sparcc[network_pred_sparcc < threshold_spacc] <- 0
 
diag(network_pred_sparcc) <- 0 #no needed

dist_matrix_spacc <- 1-network_pred_sparcc 

mantel_result_Spaac <- vegan::mantel(dist_matrix1a, dist_matrix_spacc, method = "pearson")


# only use elements in lower triangular matrix

#old code from Takemoto
real <- network_real[lower.tri(network_real)] # Aij
pred_pea <- network_pred_pea[lower.tri(network_pred_pea)] # Pearson correlation
pred_ppea <- network_pred_ppea[lower.tri(network_pred_ppea)] # Pearson partial correlation
pred_sparcc <- network_pred_sparcc[lower.tri(network_pred_sparcc)] # SparCC
str(pred_sparcc)
pred_spiec <- network_pred_spiec[lower.tri(network_pred_spiec)] #SPIEC-EASI


## Area under the Precision-Recall Curve (AUPR value)
pr.curve(pred_pea[real == 1],pred_pea[real == 0])$auc.integral # Pearson correlation
pr.curve(pred_ppea[real == 1],pred_ppea[real == 0])$auc.integral # Pearson partial correlation
pr.curve(pred_sparcc[real == 1],pred_sparcc[real == 0])$auc.integral # SparCC
pr.curve(pred_spiec[real == 1],pred_spiec[real == 0])$auc.integral #SPIEC-EASI
