
<!-- README.md is generated from README.Rmd. Please edit that file -->

# compareRhythms <img src='man/figures/logo.png' align="right" height="138.5" />

<!-- badges: start -->

[![](https://img.shields.io/badge/devel%20version-1.0.1-green.svg)](https://github.com/compareRhythms)
[![License:
GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://cran.r-project.org/web/licenses/GPL-3)
[![](https://codecov.io/gh/bharathananth/compareRhythms/branch/master/graph/badge.svg)](https://codecov.io/gh/bharathananth/compareRhythms)
[![](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![](https://img.shields.io/badge/doi-10.1111/febs.16095-yellow.svg)](https://doi.org/10.1111/febs.16095)
[![DOI](https://zenodo.org/badge/314899899.svg)](https://zenodo.org/badge/latestdoi/314899899)
<!-- badges: end -->

The goal of *compareRhythms* is to find features with altered circadian
rhythm parameters (*amplitude* and *phase*) between the control and
experimental groups.

## Installation

You can install the current version of *compareRhythms* from GitHub
with:

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.15")   # This installs bioconductor
BiocManager::install(c("SummarizedExperiment", "DESeq2", "edgeR", "limma", "rain")) # Packages needed by compareRhythms
install.packages("devtools")    # if it is not already installed
devtools::install_github("cran/npsm")   # Package archived by CRAN
devtools::install_github("cran/DODR")   # Package archived by CRAN
devtools::install_github("bharathananth/compareRhythms", build_vignettes = TRUE, dependencies = TRUE)
```

## Usage summary

The analysis is run using a single function `compareRhythms()`. To
execute this function, the three necessary ingredients are the
timeseries data, the experimental design and parameters to choose and
tune the method. The output of the function is a *data.frame* with the
IDs of the differentially rhythmic features, the category it is
classified under and optionally the rhythm parameters of the feature in
the two groups. The differential rhythmicity categories are **gain** of,
**loss** of, **change** of, or **same** rhythms (with respect to the
reference/control group).

For complete examples, please refer to the vignette titled **User
Guide** by running

``` r
library(compareRhythms)
vignette("UserGuide", package="compareRhythms")
```

## Citation

Please cite *compareRhythms* in publications as:

- Software v1.0.0 <https://doi.org/10.5281/zenodo.7699722>
- Venn diagram analysis overestimates the extent of circadian rhythm
  reprogramming. FEBS J, 289: 6605-6621.
  <https://doi.org/10.1111/febs.16095>

The R code to perform all the analyses in this publication (using
*compareRhythms*) can be found in the repository
<https://github.com/bharathananth/FEBSJ-VDA-overestimates>

## Getting help

If you encounter a bug, please file a minimal reproducible example on
[github](https://github.com/bharathananth/compareRhythms/issues).
