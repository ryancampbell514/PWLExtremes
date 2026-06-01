# Piecewise-linear geometric modelling of multivariate extremes

Code associated with <a href="https://arxiv.org/abs/2412.05195" target="_blank">the preprint</a>, by [Ryan Campbell](https://www.ryanstats.com/) and [Jennifer Wadsworth](https://www.lancaster.ac.uk/~wadswojl/#).

It is best to first familiarise yourself with the code in example_code/. The remaining scripts are for reproducing the results in the manuscript. 

## Dependencies

Please install the geometricMVE package and other dependencies using

``` r
# install geometricMVE locally:
library(remotes)
remotes::install_github("jennywadsworth/geometricMVE")

# install other required packages:
install.packages(c("geometry","evd","mvtnorm","rgl","lattice"))
```

# source in some useful functions

There are some functions requries to reproduce results in the paper, but are not yet available in the geometricMVE package. To use them, be sure to source the functions in extra-functions.R

``` r
source("extra-functions.R")
```

