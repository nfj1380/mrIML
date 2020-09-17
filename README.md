# mrIML
<img src = "https://user-images.githubusercontent.com/33707823/88988817-531ad080-d31d-11ea-8d76-f1ad0506e405.png" width="200" height="200"/>

Multivariate (multi-response) interpretable machine learning.

This package aims to enable users to build and interpret multivariate machine learning models harnessing the tidyverse (tidy model syntax in particular). This package builds off ideas from Gradient Forests (Ellis et al 2012), ecological genomic approaches (Fitzpatrick and Keller, 2014) and multi-response stacking algorithms (Xing et al 2019).

This package can be of use for any multi-response machine learning problem, but was designed to handle data common to community ecology (site by species data) and ecological genomics (individual or population by SNP loci).

# Package vignettes
Vignette for using MrIML for landscape genetic studies  [here](mrIML/Vignette_LandscapeGenetics.html)

## References
Xing, L, Lesperance, ML and Zhang, X (2020). Simultaneous prediction of multiple outcomes using revised stacking algorithms. Bioinformatics, 36, 65-72.

Fitzpatrick, M.C. & Keller, S.R. (2015) Ecological genomics meets community-level modelling of biodiversity: mapping the genomic landscape of current and future environmental adaptation. Ecology Letters 18, 1–16.

Ellis, N., Smith, S.J. and Pitcher, C.R. (2012), Gradient forests: calculating importance gradients on physical predictors. Ecology, 93: 156-168. doi:10.1890/11-0252.1
