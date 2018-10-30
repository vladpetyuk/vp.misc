vp.misc
======

## Installation
```r
source("http://bioconductor.org/biocLite.R")
options(repos=biocinstallRepos(character()))
install.packages("devtools")
library("devtools")
install_github("vladpetyuk/vp.misc")
```

## Description
This is a collection of misc custom functions.  They are primarily geared towards proteomics/perhaps-other-omics data analysis using `MSnSet` or `ExpressionSet` as data containers.  
* heatmaps
* colorschemes
* nested linear models
* wrappers for the machine learning routines
* custom imputation
* reading of MaxQuant and SkyLine data
* handling PTM mapping
