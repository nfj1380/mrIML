#####################################################################
#PART 1: Microbiome simulation
#
# Contact:
# Kazuhiro Takemoto (takemoto@bio.kyutech.ac.jp)
#####################################################################

library(seqtime) # https://github.com/hallucigenia-sparsa/seqtime
library(SpiecEasi) # https://github.com/zdk123/SpiecEasi
library(igraph) # https://igraph.org/r/
library(ppcor) # https://cran.r-project.org/web/packages/ppcor/index.html
library(PRROC) # https://cran.r-project.org/package=PRROC
# load generateM_specific_type function
source("generateM_specific_type.R")

# library(devtools)  
# install_github("hallucigenia-sparsa/seqtime")  
# library(seqtime) 

nn <- 40 # network size (n)
k_ave <- 2 # average degree (<k>). Could base this off real world nets?

## Genrate an interaction matrix (Mij)
obj <- generateM_specific_type(nn,k_ave,type.network="random",type.interact="mix",interact.str.max=0.8,mix.compt.ratio=0.5)
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
network_real <- obj[[1]]
# interaction matrix (Mij) for the GLV model
M <- obj[[2]]
Mdf <-as.data.frame( obj[[2]])
colnames(Mdf) <- row.names(Mdf)
#convert to igraph
diag(M) <- 0
diag(Mdf) <- 0
g <- graph.adjacency(M, mode = "undirected", weighted = TRUE)#matching Y data

# Set edge colors based on positive (red) and negative (blue) weights
edge_colors <- ifelse(E(g)$weight > 0, "red", "blue")

# Set edge widths based on the absolute values of edge weights
edge_widths <- abs(E(g)$weight) * 5  # Adjust the scaling factor as needed

plot(
  g,
  #vertex.label = NA,
  edge.color = edge_colors,
  edge.width = edge_widths
)
gg <- ggnetwork(g)

# Plot the graph
ggplot(gg, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(aes(color = edge_colors, linewidth = (edge_widths)), curvature = 0.2,
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


## plot population dynamics
y <- rpois(nn,lambda=100)
r <- runif(nn)
res <- glv(nn, M, r, y)
tsplot(10*res[,20:1000],time.given =T)

## Generate dataset on species abundance using the GLV model. 300 samples
data <- as.data.frame(generateDataSet(300, M, count = nn*100, mode = 4)) #mode follows from Hirano and Takemoto 2019
# relative abundance
data_relative <- as.data.frame(t(data) / apply(data,2,sum))
str(data_relative)
 data_mean <- colMeans(data_relative )

# Transform the data to meet the condition
data_pa  <- data_relative %>%
  mutate(across(everything(), ~ ifelse(. < data_mean, 0, 1)))# %>%
  #mutate_all(factor)


##############################################################
#MrIML analysis
##############################################################

Y <- filterRareCommon(data_pa, lower=0.1, higher=0.9)  #9 lost

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
                              prop=0.5, k=5, racing=F) #racing was having issues

ModelPerf <- mrIMLperformance(yhats_rf_sim , Model=model_rf, Y=Y, mode='classification')
m <- ModelPerf[[1]]

bs_sim <- mrBootstrap(yhats=yhats_rf_sim, Y=Y, num_bootstrap = 20, alpha = 0.05, ice=F)

bs_impVIa <- mrVI_bootstrap(mrBootstrap_obj=bs_sim, ModelPerf=ModelPerf, 
                            threshold=0.9,  X=X1, Y=Y, global_top_var=10,
                            local_top_var=5)
vi_obj <- bs_impVIa[[1]]#data for posterity
bs_impVIa[[2]]#combined plot

#create bootstrapped pdps. Global_top5 allows just to plot the # most important features

pds <- mrPD_bootstrap(mrBootstrap_obj=bs_sim, vi_obj=bs_impVIa, X=X1, Y=Y,
                      target='sp32', global_top_var=5)
pds[[1]] #data
pds[[2]]#plot

############################################
#Network
##############################################
assoc_net<- mrCoOccurNet_bootstrap (mrPD_obj=pds , Y=Y)

assoc_net_filtered <-  assoc_net %>% 
  filter(mean_strength > 0.025) #seems like the sweetspot? worse a 0.1

#to make everthing comparable
transformed_strength <- assoc_net_filtered  %>%
  mutate(mean_strength_dir = ifelse(direction == "negative", -mean_strength, mean_strength))


#convert to igraph
g <- graph_from_data_frame(assoc_net_filtered, directed=FALSE, vertices=names(Y)) #matching Y data

#need to fix this - reformat the dataframe
#plot(g)

E(g)$Value <- assoc_net_filtered$mean_strength###chnge this as needed
E(g)$Color <- ifelse(assoc_net_filtered$direction == "negative", "blue", "red")
E(g)$Value2 <- transformed_strength$mean_strength_dir #
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


###################################################
#How well does MrIML do?
###################################################

mrIML_matrix <- get.adjacency(g, attr = "Value2", sparse = FALSE)
mrIML_matrix

#need to make sure everything matches for each network

# Remove 'sp' prefix from row names
new_row_names <- sub("^sp", "", rownames(mrIML_matrix))

# Set the new row names
rownames(mrIML_matrix) <- new_row_names

matching_indices <- match(rownames(mrIML_matrix), rownames(Mdf))
adj_M_reordered <- M[matching_indices, matching_indices]

dist_matrix1 <- 1 - adj_M_reordered
dist_matrix2 <- 1 - mrIML_matrix

#mantel test to see if these adjacency matrices are different or not
mantel_result <- vegan::mantel(dist_matrix1, dist_matrix2, method = "pearson")
mantel_result #not bad for a first attempt 0.52 and significant

###################################################
## Inffering ecological associations (example)
# based on Pearson correlation
network_pred_pea <- abs(cor(t(data_relative)))
# based on Pearson partial correlation
network_pred_ppea <- abs(pcor(t(data_relative))$estimate)
# based on SparCC
# Note that in this study we used SparCC python module,
# not the follwing R wrapper function (see Methods section in the main text).
# (But, similar results are obtained)
network_pred_sparcc <- abs(sparcc(t(data))$Cor)
# based on SPIEC-EASI
network_pred_spiec <- spiec.easi(t(data),method='mb')
network_pred_spiec <- as.matrix(getOptMerge(network_pred_spiec))


## Evaluating evaluating co-occurrence network performance
# only use elements in lower triangular matrix
real <- network_real[lower.tri(network_real)] # Aij
pred_pea <- network_pred_pea[lower.tri(network_pred_pea)] # Pearson correlation
pred_ppea <- network_pred_ppea[lower.tri(network_pred_ppea)] # Pearson partial correlation
pred_sparcc <- network_pred_sparcc[lower.tri(network_pred_sparcc)] # SparCC
pred_spiec <- network_pred_spiec[lower.tri(network_pred_spiec)] #SPIEC-EASI


## Area under the Precision-Recall Curve (AUPR value)
pr.curve(pred_pea[real == 1],pred_pea[real == 0])$auc.integral # Pearson correlation
pr.curve(pred_ppea[real == 1],pred_ppea[real == 0])$auc.integral # Pearson partial correlation
pr.curve(pred_sparcc[real == 1],pred_sparcc[real == 0])$auc.integral # SparCC
pr.curve(pred_spiec[real == 1],pred_spiec[real == 0])$auc.integral #SPIEC-EASI
