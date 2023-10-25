
# Overview

`mcbrnet` is a collection of functions to perform ecological simulations
in spatial networks. Main functions in this package can: (1) generate a
random branching network with the specified number of patches and
probability of branching (`brnet()`); and (2) simulate
metacommunity/foodweb dynamics in space (`mcsim()`, `igpsim()`,
`sglv()`). See **Articles** for instruction.

# Citation

Current best citation for this package is:

[Terui A, Kim S, Dolph CL, Kadoya T, Miyazaki Y. (2021) Emergent dual
scaling of riverine biodiversity. *Proceedings of the National Academy
of Sciences* 118: e2105574118](https://doi.org/10.1073/pnas.2105574118)

# Installation

The `mcbrnet` package can be installed with the following script:

``` r
#install.packages("remotes")
remotes::install_github("aterui/mcbrnet")
library(mcbrnet)
```

# Funding

This material is based upon work supported by the National Science
Foundation through the Division of Environmental Biology (DEB 2015634).
