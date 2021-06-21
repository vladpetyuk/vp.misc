The package moved to [PNNL-Comp-Mass-Spec](https://github.com/PNNL-Comp-Mass-Spec) as [MSnSet.utils](https://github.com/PNNL-Comp-Mass-Spec/MSnSet.utils). This repository may or may not be maintained.



vp.misc
======

## Installation
Make sure package `devtools` is installed.

To install `vp.misc` run
```r
devtools::install_github("vladpetyuk/vp.misc")
```
By default `install_github` does not compile vignettes. Thus to compile with vignettes use the following options.
```r
devtools::install_github("vladpetyuk/vp.misc", build = TRUE, build_opts = c("--no-resave-data", "--no-manual"))
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
