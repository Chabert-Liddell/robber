
<!-- README.md is generated from README.Rmd. Please edit that file -->

# robber

<!-- badges: start -->
[![](https://www.r-pkg.org/badges/version/robber?color=green)](https://cran.r-project.org/package=robber)
[![R build
status](https://github.com/Chabert-Liddell/robber/workflows/R-CMD-check/badge.svg)](https://github.com/Chabert-Liddell/robber/actions)
[![Codecov test
coverage](https://codecov.io/gh/Chabert-Liddell/robber/branch/master/graph/badge.svg)](https://app.codecov.io/gh/Chabert-Liddell/robber?branch=master)
<!-- badges: end -->

ROBustness using Blockmodel for Ecological network R package

Implementation of a variety of methods to compute
    the robustness of ecological interaction networks with binary interactions 
    as described in [Chabert-Liddell, Barbillon and Donnet (2022)](https://doi.org/10.1002/env.2709). In particular, using the Stochastic 
    Block Model and its bipartite counterpart, the Latent Block Model to put a 
    parametric model on the network, allows the comparison of the robustness of 
    networks differing in species richness and number of interactions. It also
    deals with networks that are partially sampled and/or with missing values. 




## Installation

You can install the latest cran release of robber with:
``` r
install.packages("robber")
```

or the development version of robber
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Chabert-Liddell/robber")
```
## References

Chabert-Liddell, S.-C., Barbillon, P., & Donnet, S. (2022). Impact of the mesoscale structure of a bipartite ecological interaction network on its robustness through a probabilistic modeling. Environmetrics, 33( 2), e2709. https://doi.org/10.1002/env.2709
