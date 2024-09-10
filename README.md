mcbrnet
================

# Overview

The R package `mcbrnet` provides functions to perform ecological
simulations in spatial networks. Main functions in this package can: (1)
generate a random branching network with the specified number of patches
and probability of branching (`brnet()`); and (2) simulate
metacommunity/foodweb dynamics in space (`mcsim()`, `igpsim()`). See
**Articles** for instruction.

# Citation

Current best citations for this package is:

[Terui A, Kim S, Dolph CL, Kadoya T, Miyazaki Y. (2021) Emergent dual
scaling of riverine biodiversity. *Proceedings of the National Academy
of Sciences* 118: e2105574118](https://doi.org/10.1073/pnas.2105574118)

[Pomeranz JPF, Finlay JF, Terui A. (2023) Ecosystem size and complexity
as extrinsic drivers of food chain length in branching ecosystems.
*Ecosphere* 14(8):
e4648](https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecs2.4648)

# Installation

The latest stable version of `mcbrnet` package can be installed with:

``` r
#install.packages("remotes")
remotes::install_github("aterui/mcbrnet@1.4.1")
library(mcbrnet)
```

The development version of `mcbrnet` package can be installed with:

``` r
#install.packages("remotes")
remotes::install_github("aterui/mcbrnet")
library(mcbrnet)
```

# Funding

This material is based upon work supported by the National Science
Foundation through the Division of Environmental Biology (DEB 2015634).
