
# Overview

`mcbrnet` is a collection of functions to perform metacommunity
simulation in branching networks. Main functions include:

- `brnet()` generates a random branching network with the specified
  number of patches and probability of branching. The function returns
  adjacency and distance matrices, hypothetical
  environmental/disturbance values at each patch, and the number of
  patches upstream (akin to the watershed area in river networks).

- `mcsim()` simulates metacommunity dynamics. This function follows a
  general framework proposed by [Thompson et
  al. (2020)](https://doi.org/10.1111/ele.13568). However, it has
  several unique features that are detailed in [Terui et al
  (2021)](https://doi.org/10.1073/pnas.2105574118).

- `igpsim()` simulates three-species meta-food web dynamics with
  intraguild predation. This function shares many features with
  `mcsim()`. Currently under development.

- `ggbrnet()` is a wrapper of `ggraph` functions for easy visualization
  of a network produced by `brnet()` .

- `ptsource()` is a function to simulate propagation of environmental
  pollutants in a river network.

- `frgm()` is a function to simulate fragmentation in a river network.

- `adjtodist()` converts an adjacency matrix to a distance matrix.

See Articles for instruction.

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
