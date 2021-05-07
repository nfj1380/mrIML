
<!-- README.md is generated from README.Rmd. Please edit that file -->
# mrIML: Multivariate (multi-response) interpretable machine learning <a href='https://nfj1380.github.io/mrIML/index.html'><img src="man/figures/logo.png" align="right" height="139"/></a>

<!-- badges: start -->
![GitHub R package version](https://img.shields.io/github/r-package/v/nfj1380/mrIML?logo=github&logoColor=%2300ff37&style=flat-square) ![GitHub contributors](https://img.shields.io/github/contributors/nfj1380/mrIML?style=flat-square) ![GitHub all releases](https://img.shields.io/github/downloads/nfj1380/mriml/total?style=flat-square) ![GitHub last commit](https://img.shields.io/github/last-commit/nfj1380/mrIML?style=flat-square) [![R-CMD-check](https://github.com/nfj1380/mrIML/workflows/R-CMD-check/badge.svg)](https://github.com/nfj1380/mrIML/actions) <!-- badges: end -->

This package aims to enable users to build and interpret multivariate machine learning models harnessing the tidyverse (tidy model syntax in particular). This package builds off ideas from Gradient Forests (Ellis et al 2012), ecological genomic approaches (Fitzpatrick and Keller, 2014) and multi-response stacking algorithms (Xing et al 2019).

This package can be of use for any multi-response machine learning problem, but was designed to handle data common to community ecology (site by species data) and ecological genomics (individual or population by SNP loci).

## Installation

Install the stable version of the package:

``` r
#install.packages("devtools")
devtools::install_github('nfj1380/mrIML')
```

## Quick start

**mrIML** is designed to be used with a single function call or to be used in an ad-hoc fashion via individual function calls. In the following section we give an overview of the simple use case. For more on using each function see the [function documentation](xx). The core functions for both regression and classification are: [`mrIMLpredicts`](xx), [`mrIMLperformance`](xx), and [`mrInteractions`](xx),for plotting and visualization [`mrVip`](xx), [`mrFlashlight`](xx), and[`plot_vi`](xx). Estimating the interactions alone can be substantially computationally demanding depending on the number of outcomes you want to test. The first step to using the package is to load it as follows.

``` r
library(mrIML)
#other package needed:
library(vip); library(tidymodels); library(randomForest);  library(caret); library(gbm);
#> Warning: package 'vip' was built under R version 4.0.5
#> Warning: package 'randomForest' was built under R version 4.0.5
#> Warning: package 'caret' was built under R version 4.0.5
#> Warning: package 'gbm' was built under R version 4.0.5
library(tidyverse);library(parallel); library(doParallel); library(themis); library(viridis);
#> Warning: package 'doParallel' was built under R version 4.0.5
#> Warning: package 'themis' was built under R version 4.0.5
library(janitor); library(hrbrthemes); library(xgboost); library(vegan);library(flashlight);
#> Warning: package 'janitor' was built under R version 4.0.5
#> Warning: package 'hrbrthemes' was built under R version 4.0.5
#> Warning: package 'xgboost' was built under R version 4.0.5
#> Warning: package 'vegan' was built under R version 4.0.5
#> Warning: package 'permute' was built under R version 4.0.5
#> Warning: package 'flashlight' was built under R version 4.0.5
library(ggrepel); library(parsnip);library(rsample); library(workflows)
#> Warning: package 'ggrepel' was built under R version 4.0.5
```

## Model component

Now all the data is loaded and ready to go we can formulate the model using tidymodel syntax. In this case we have binary data (SNP presence/absence at each loci) but the data could also be counts or continuous (the set\_model argument would be "regression" instead of "classification"). The user can specify any model from the tidymodel universe as 'model 1' (see <https://www.tidymodels.org/find/> for details). However, we have done most of our testing on random forests (rf), xgr boost and glms (generalized linear models). Here we will specify a random forest classification model as the model applied to each response.

``` r
model1 <- 
  rand_forest(trees = 100, mode = "classification") %>% #100 trees are set for brevity
  set_engine("ranger", importance = c("impurity","impurity_corrected")) %>%# select the engine/package that underlies the model
  set_mode("classification")# choose either the continuous "regression" or binary "classification" mode
```

### \[mrIMLpredicts\]

This function represents the core functionality of the package and includes results reporting, plotting and optional saving. It requires a data frame of X t( the snp data for example) and Y represented by the covariates or features.

Load example data (cite) data from `{mrIML}`.

``` r
data <- gfData[1:20]
head(data)
#>         x     y bio_1 bio_10 bio_15 bio_18 bio_2 bio_7 elevation     MEM.1      MEM.2      MEM.3
#> 1 -101.57 49.14    31    182     63    189   131   489       498  24.55053   7.444415  -9.018635
#> 2 -147.89 69.10  -107     85     60     89   110   507       400 -14.18000 -22.945314 -11.845836
#> 3 -148.51 63.39   -52     90     63    196   102   399      1040 -15.77654 -27.420568 -13.038813
#> 4  -98.35 54.51   -12    159     52    204   103   510       213  20.73482   8.715802 -19.187171
#> 5 -112.79 53.33    18    150     66    234   119   427       778  18.51313  -3.488470  24.300653
#> 6 -146.17 64.59   -43    135     65    154   109   493       460 -17.18802 -31.070507 -13.633361
#>         MEM.4 REFERENCE_X80928 REFERENCE_X172196 REFERENCE_X172210 REFERENCE_X172625
#> 1  0.06879434        0.4000000        0.03333333         0.2692308         0.4230769
#> 2 21.59719769        0.3611111        0.00000000         0.2222222         0.4705882
#> 3 14.65433486        0.4333333        0.00000000         0.1666667         0.1666667
#> 4 -5.13771830        0.3846154        0.07142857         0.4230769         0.4615385
#> 5  7.73910268        0.4285714        0.00000000         0.4000000         0.4000000
#> 6  5.78938201        0.3333333        0.00000000         0.1785714         0.3846154
#>   REFERENCE_X177996 REFERENCE_X180626 REFERENCE_X181299
#> 1        0.00000000        0.03846154        0.06666667
#> 2        0.00000000        0.05555556        0.00000000
#> 3        0.03571429        0.16666667        0.00000000
#> 4        0.10714286        0.00000000        0.00000000
#> 5        0.06666667        0.19230769        0.06666667
#> 6        0.00000000        0.08333333        0.00000000
```

``` r
# Define set of features
FeaturesnoNA<-Features[complete.cases(Features), ]
Y <- FeaturesnoNA #for simplicity
# Define set the outcomes of interst
fData <- filterRareCommon (Responsedata, lower=0.4, higher=0.7) 
X <- fData #

yhats <- mrIMLpredicts(X=X,Y=Y, model1=model1, balance_data='no', mod='classification', parallel = TRUE)
#save(yhats, file='logreg_model')
ModelPerf <- mrIMLperformance(yhats, model1, X=X) #
ModelPerf[[2]] #this measures performance across all loci.
#> [1] 0.5517241
```

## Plotting

``` r
VI <- mrVip(yhats, Y=Y) 
plot_vi(VI=VI,  X=X,Y=Y, modelPerf=ModelPerf, cutoff= 0, plot.pca='yes') #the cutoff reduces the number of individual models printed in the second plot. 
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" />

    #> Press [enter] to plot individual variable importance summaries

<img src="man/figures/README-unnamed-chunk-7-2.png" width="100%" />

    #> Press [enter] to plot the importance PCA plot

<img src="man/figures/README-unnamed-chunk-7-3.png" width="100%" />

## Effect of a feature on genetic change

We also wrap some flashlight functionality to visualize the marginal (i.e. partial dependencies) or conditional (accumulated local effects) effect of a feature on genetic change. Partial dependencies take longer to calculate and are more sensitive to correlated features

``` r
flashlightObj <- mrFlashlight(yhats, X, Y, response = "multi", model='classification')

#plot prediction scatter for all responses. Gets busy with 
plot(light_scatter(flashlightObj, v = "Forest", type = "predicted"))
```

<img src="man/figures/README-unnamed-chunk-9-1.png" width="100%" />

``` r

#plots everything on one plot (partial dependency, ALE, scatter)
plot(light_effects(flashlightObj, v = "Grassland"), use = "all")
```

<img src="man/figures/README-unnamed-chunk-9-2.png" width="100%" />

``` r


#profileData_pd <- light_profile(flashlightObj,  v = "Grassland")

#mrProfileplot(profileData_pd , sdthresh =0.05) #sdthresh removes responses from the first plot that do not vary with the feature

profileData_ale <- light_profile(flashlightObj, v = "Grassland", type = "ale") #acumulated local effects

mrProfileplot(profileData_ale , sdthresh =0.01)
```

<img src="man/figures/README-unnamed-chunk-9-3.png" width="100%" /><img src="man/figures/README-unnamed-chunk-9-4.png" width="100%" />

``` r
#the second plot is the cumulative turnover function
```

## Interacting predictors or features

Finally, we can assess how features interact overall to shape genetic change. Be warned this is memory intensive. Future updates to this package will enable users to visualize these interactions and explore them in more detail using 2D ALE plots for example.

``` r

#interactions <-mrInteractions(yhats, X, Y,  mod='classification') #this is computationally intensive so multicores are needed. If stopped prematurely - have to reload things

#mrPlot_interactions(interactions, X,Y, top_ranking = 2, top_response=2)
```

## References

Xing, L, Lesperance, ML and Zhang, X (2020). Simultaneous prediction of multiple outcomes using revised stacking algorithms. Bioinformatics, 36, 65-72.

Fitzpatrick, M.C. & Keller, S.R. (2015) Ecological genomics meets community-level modelling of biodiversity: mapping the genomic landscape of current and future environmental adaptation. Ecology Letters 18, 1–16.

Ellis, N., Smith, S.J. and Pitcher, C.R. (2012), Gradient forests: calculating importance gradients on physical predictors. Ecology, 93: 156-168. <doi:10.1890/11-0252.1>

## Session info

``` r
devtools::session_info()
#> - Session info -----------------------------------------------------------------------------------
#>  setting  value                       
#>  version  R version 4.0.4 (2021-02-15)
#>  os       Windows 10 x64              
#>  system   x86_64, mingw32             
#>  ui       RTerm                       
#>  language (EN)                        
#>  collate  English_Australia.1252      
#>  ctype    English_Australia.1252      
#>  tz       Australia/Hobart            
#>  date     2021-05-07                  
#> 
#> - Packages ---------------------------------------------------------------------------------------
#>  package         * version    date       lib source                        
#>  assertthat        0.2.1      2019-03-21 [1] CRAN (R 4.0.3)                
#>  backports         1.2.1      2020-12-09 [1] CRAN (R 4.0.3)                
#>  BBmisc            1.11       2017-03-10 [1] CRAN (R 4.0.5)                
#>  broom           * 0.7.5      2021-02-19 [1] CRAN (R 4.0.4)                
#>  cachem            1.0.4      2021-02-13 [1] CRAN (R 4.0.4)                
#>  callr             3.5.1      2020-10-13 [1] CRAN (R 4.0.3)                
#>  caret           * 6.0-86     2020-03-20 [1] CRAN (R 4.0.5)                
#>  cellranger        1.1.0      2016-07-27 [1] CRAN (R 4.0.3)                
#>  checkmate         2.0.0      2020-02-06 [1] CRAN (R 4.0.3)                
#>  class             7.3-18     2021-01-24 [2] CRAN (R 4.0.4)                
#>  cli               2.3.1      2021-02-23 [1] CRAN (R 4.0.4)                
#>  cluster           2.1.0      2019-06-19 [2] CRAN (R 4.0.4)                
#>  codetools         0.2-18     2020-11-04 [2] CRAN (R 4.0.4)                
#>  colorspace        2.0-0      2020-11-11 [1] CRAN (R 4.0.3)                
#>  cowplot           1.1.1      2020-12-30 [1] CRAN (R 4.0.4)                
#>  crayon            1.4.1      2021-02-08 [1] CRAN (R 4.0.3)                
#>  data.table        1.14.0     2021-02-21 [1] CRAN (R 4.0.4)                
#>  DBI               1.1.1      2021-01-15 [1] CRAN (R 4.0.3)                
#>  dbplyr            2.1.0      2021-02-03 [1] CRAN (R 4.0.3)                
#>  desc              1.2.0      2018-05-01 [1] CRAN (R 4.0.3)                
#>  devtools          2.3.2      2020-09-18 [1] CRAN (R 4.0.3)                
#>  dials           * 0.0.9      2020-09-16 [1] CRAN (R 4.0.4)                
#>  DiceDesign        1.9        2021-02-13 [1] CRAN (R 4.0.4)                
#>  digest            0.6.27     2020-10-24 [1] CRAN (R 4.0.3)                
#>  doParallel      * 1.0.16     2020-10-16 [1] CRAN (R 4.0.5)                
#>  dplyr           * 1.0.4      2021-02-02 [1] CRAN (R 4.0.3)                
#>  ellipsis          0.3.1      2020-05-15 [1] CRAN (R 4.0.3)                
#>  evaluate          0.14       2019-05-28 [1] CRAN (R 4.0.3)                
#>  extrafont         0.17       2014-12-08 [1] CRAN (R 4.0.3)                
#>  extrafontdb       1.0        2012-06-11 [1] CRAN (R 4.0.3)                
#>  fansi             0.4.2      2021-01-15 [1] CRAN (R 4.0.3)                
#>  farver            2.0.3      2020-01-16 [1] CRAN (R 4.0.3)                
#>  fastmap           1.1.0      2021-01-25 [1] CRAN (R 4.0.3)                
#>  fastmatch         1.1-0      2017-01-28 [1] CRAN (R 4.0.3)                
#>  flashlight      * 0.7.5      2021-02-13 [1] CRAN (R 4.0.5)                
#>  FNN               1.1.3      2019-02-15 [1] CRAN (R 4.0.4)                
#>  forcats         * 0.5.1      2021-01-27 [1] CRAN (R 4.0.3)                
#>  foreach         * 1.5.1      2020-10-15 [1] CRAN (R 4.0.3)                
#>  fs                1.5.0      2020-07-31 [1] CRAN (R 4.0.3)                
#>  furrr             0.2.2      2021-01-29 [1] CRAN (R 4.0.4)                
#>  future            1.21.0     2020-12-10 [1] CRAN (R 4.0.3)                
#>  gbm             * 2.1.8      2020-07-15 [1] CRAN (R 4.0.5)                
#>  gdtools           0.2.3      2021-01-06 [1] CRAN (R 4.0.5)                
#>  generics          0.1.0      2020-10-31 [1] CRAN (R 4.0.3)                
#>  ggplot2         * 3.3.3      2020-12-30 [1] CRAN (R 4.0.3)                
#>  ggrepel         * 0.9.1      2021-01-15 [1] CRAN (R 4.0.5)                
#>  globals           0.14.0     2020-11-22 [1] CRAN (R 4.0.3)                
#>  glue              1.4.2      2020-08-27 [1] CRAN (R 4.0.3)                
#>  gower             0.2.2      2020-06-23 [1] CRAN (R 4.0.3)                
#>  GPfit             1.0-8      2019-02-08 [1] CRAN (R 4.0.4)                
#>  gridExtra         2.3        2017-09-09 [1] CRAN (R 4.0.3)                
#>  gtable            0.3.0      2019-03-25 [1] CRAN (R 4.0.3)                
#>  hardhat           0.1.5      2020-11-09 [1] CRAN (R 4.0.4)                
#>  haven             2.3.1      2020-06-01 [1] CRAN (R 4.0.3)                
#>  highr             0.8        2019-03-20 [1] CRAN (R 4.0.3)                
#>  hms               1.0.0      2021-01-13 [1] CRAN (R 4.0.3)                
#>  hrbrthemes      * 0.8.0      2020-03-06 [1] CRAN (R 4.0.5)                
#>  htmltools         0.5.1.1    2021-01-22 [1] CRAN (R 4.0.3)                
#>  httr              1.4.2      2020-07-20 [1] CRAN (R 4.0.3)                
#>  infer           * 0.5.4      2021-01-13 [1] CRAN (R 4.0.4)                
#>  ipred             0.9-11     2021-03-12 [1] CRAN (R 4.0.4)                
#>  iterators       * 1.0.13     2020-10-15 [1] CRAN (R 4.0.3)                
#>  janitor         * 2.1.0      2021-01-05 [1] CRAN (R 4.0.5)                
#>  jsonlite          1.7.2      2020-12-09 [1] CRAN (R 4.0.3)                
#>  knitr             1.31       2021-01-27 [1] CRAN (R 4.0.3)                
#>  labeling          0.4.2      2020-10-20 [1] CRAN (R 4.0.3)                
#>  lattice         * 0.20-41    2020-04-02 [2] CRAN (R 4.0.4)                
#>  lava              1.6.9      2021-03-11 [1] CRAN (R 4.0.4)                
#>  lhs               1.1.1      2020-10-05 [1] CRAN (R 4.0.4)                
#>  lifecycle         1.0.0      2021-02-15 [1] CRAN (R 4.0.4)                
#>  listenv           0.8.0      2019-12-05 [1] CRAN (R 4.0.3)                
#>  lubridate         1.7.9.2    2020-11-13 [1] CRAN (R 4.0.3)                
#>  magrittr          2.0.1      2020-11-17 [1] CRAN (R 4.0.3)                
#>  MASS              7.3-53     2020-09-09 [2] CRAN (R 4.0.4)                
#>  Matrix            1.3-2      2021-01-06 [2] CRAN (R 4.0.4)                
#>  memoise           2.0.0      2021-01-26 [1] CRAN (R 4.0.3)                
#>  MetricsWeighted   0.5.1      2020-04-18 [1] CRAN (R 4.0.5)                
#>  mgcv              1.8-33     2020-08-27 [2] CRAN (R 4.0.4)                
#>  mlr               2.19.0     2021-02-22 [1] CRAN (R 4.0.5)                
#>  modeldata       * 0.1.0      2020-10-22 [1] CRAN (R 4.0.4)                
#>  ModelMetrics      1.2.2.2    2020-03-17 [1] CRAN (R 4.0.5)                
#>  modelr            0.1.8      2020-05-19 [1] CRAN (R 4.0.3)                
#>  mrIML           * 0.1        2021-04-13 [1] Github (nfj1380/mrIML@57678b8)
#>  munsell           0.5.0      2018-06-12 [1] CRAN (R 4.0.3)                
#>  nlme              3.1-152    2021-02-04 [2] CRAN (R 4.0.4)                
#>  nnet              7.3-15     2021-01-24 [2] CRAN (R 4.0.4)                
#>  parallelly        1.23.0     2021-01-04 [1] CRAN (R 4.0.3)                
#>  parallelMap       1.5.0      2020-03-26 [1] CRAN (R 4.0.5)                
#>  ParamHelpers      1.14       2020-03-24 [1] CRAN (R 4.0.5)                
#>  parsnip         * 0.1.5      2021-01-19 [1] CRAN (R 4.0.4)                
#>  permute         * 0.9-5      2019-03-12 [1] CRAN (R 4.0.5)                
#>  pillar            1.5.1      2021-03-05 [1] CRAN (R 4.0.4)                
#>  pkgbuild          1.2.0      2020-12-15 [1] CRAN (R 4.0.3)                
#>  pkgconfig         2.0.3      2019-09-22 [1] CRAN (R 4.0.3)                
#>  pkgload           1.1.0      2020-05-29 [1] CRAN (R 4.0.3)                
#>  plyr              1.8.6      2020-03-03 [1] CRAN (R 4.0.3)                
#>  prettyunits       1.1.1      2020-01-24 [1] CRAN (R 4.0.3)                
#>  pROC              1.17.0.1   2021-01-13 [1] CRAN (R 4.0.4)                
#>  processx          3.4.5      2020-11-30 [1] CRAN (R 4.0.3)                
#>  prodlim           2019.11.13 2019-11-17 [1] CRAN (R 4.0.4)                
#>  ps                1.5.0      2020-12-05 [1] CRAN (R 4.0.3)                
#>  purrr           * 0.3.4      2020-04-17 [1] CRAN (R 4.0.3)                
#>  R6                2.5.0      2020-10-28 [1] CRAN (R 4.0.3)                
#>  randomForest    * 4.6-14     2018-03-25 [1] CRAN (R 4.0.5)                
#>  ranger            0.12.1     2020-01-10 [1] CRAN (R 4.0.4)                
#>  RANN              2.6.1      2019-01-08 [1] CRAN (R 4.0.5)                
#>  Rcpp              1.0.6      2021-01-15 [1] CRAN (R 4.0.3)                
#>  readr           * 1.4.0      2020-10-05 [1] CRAN (R 4.0.3)                
#>  readxl            1.3.1      2019-03-13 [1] CRAN (R 4.0.3)                
#>  recipes         * 0.1.15     2020-11-11 [1] CRAN (R 4.0.4)                
#>  remotes           2.2.0      2020-07-21 [1] CRAN (R 4.0.3)                
#>  reprex            1.0.0      2021-01-27 [1] CRAN (R 4.0.3)                
#>  reshape2          1.4.4      2020-04-09 [1] CRAN (R 4.0.3)                
#>  rlang             0.4.10     2020-12-30 [1] CRAN (R 4.0.4)                
#>  rmarkdown         2.7        2021-02-19 [1] CRAN (R 4.0.4)                
#>  ROSE              0.0-3      2014-07-15 [1] CRAN (R 4.0.5)                
#>  rpart             4.1-15     2019-04-12 [2] CRAN (R 4.0.4)                
#>  rpart.plot        3.0.9      2020-09-17 [1] CRAN (R 4.0.5)                
#>  rprojroot         2.0.2      2020-11-15 [1] CRAN (R 4.0.3)                
#>  rsample         * 0.0.9      2021-02-17 [1] CRAN (R 4.0.4)                
#>  rstudioapi        0.13       2020-11-12 [1] CRAN (R 4.0.3)                
#>  Rttf2pt1          1.3.8      2020-01-10 [1] CRAN (R 4.0.3)                
#>  rvest             0.3.6      2020-07-25 [1] CRAN (R 4.0.3)                
#>  scales          * 1.1.1      2020-05-11 [1] CRAN (R 4.0.3)                
#>  sessioninfo       1.1.1      2018-11-05 [1] CRAN (R 4.0.3)                
#>  snakecase         0.11.0     2019-05-25 [1] CRAN (R 4.0.5)                
#>  stringi           1.5.3      2020-09-09 [1] CRAN (R 4.0.3)                
#>  stringr         * 1.4.0      2019-02-10 [1] CRAN (R 4.0.3)                
#>  survival          3.2-7      2020-09-28 [2] CRAN (R 4.0.4)                
#>  systemfonts       1.0.1      2021-02-09 [1] CRAN (R 4.0.4)                
#>  testthat          3.0.2      2021-02-14 [1] CRAN (R 4.0.4)                
#>  themis          * 0.1.3      2020-11-12 [1] CRAN (R 4.0.5)                
#>  tibble          * 3.1.0      2021-02-25 [1] CRAN (R 4.0.4)                
#>  tidymodels      * 0.1.2      2020-11-22 [1] CRAN (R 4.0.4)                
#>  tidyr           * 1.1.2      2020-08-27 [1] CRAN (R 4.0.3)                
#>  tidyselect        1.1.0      2020-05-11 [1] CRAN (R 4.0.3)                
#>  tidyverse       * 1.3.0      2019-11-21 [1] CRAN (R 4.0.3)                
#>  timeDate          3043.102   2018-02-21 [1] CRAN (R 4.0.4)                
#>  tune            * 0.1.3      2021-02-28 [1] CRAN (R 4.0.4)                
#>  unbalanced        2.0        2015-06-26 [1] CRAN (R 4.0.5)                
#>  usethis           2.0.1      2021-02-10 [1] CRAN (R 4.0.3)                
#>  utf8              1.1.4      2018-05-24 [1] CRAN (R 4.0.3)                
#>  vctrs             0.3.6      2020-12-17 [1] CRAN (R 4.0.3)                
#>  vegan           * 2.5-7      2020-11-28 [1] CRAN (R 4.0.5)                
#>  vip             * 0.3.2      2020-12-17 [1] CRAN (R 4.0.5)                
#>  viridis         * 0.5.1      2018-03-29 [1] CRAN (R 4.0.3)                
#>  viridisLite     * 0.3.0      2018-02-01 [1] CRAN (R 4.0.3)                
#>  withr             2.4.1      2021-01-26 [1] CRAN (R 4.0.3)                
#>  workflows       * 0.2.2      2021-03-10 [1] CRAN (R 4.0.4)                
#>  xfun              0.21       2021-02-10 [1] CRAN (R 4.0.3)                
#>  xgboost         * 1.3.2.1    2021-01-18 [1] CRAN (R 4.0.5)                
#>  xml2              1.3.2      2020-04-23 [1] CRAN (R 4.0.3)                
#>  yaml              2.2.1      2020-02-01 [1] CRAN (R 4.0.3)                
#>  yardstick       * 0.0.8      2021-03-28 [1] CRAN (R 4.0.4)                
#> 
#> [1] C:/Users/Nick FJ/Documents/R/win-library/4.0
#> [2] C:/Program Files/R/R-4.0.4/library
```
