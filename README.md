mcbrnet
================

# Overview

Main functions in `mcbrnet` ;

-   `brnet`: Function `brnet()` generates a random branching network
    with the specified number of patches and probability of branching.
    The function returns adjacency and distance matrices, hypothetical
    environmental/disturbance values at each patch, and the number of
    patches upstream (akin to the watershed area in river networks). The
    output may be used in function `mcsim()` to simulate metacommunity
    dynamics in a branching network.

-   `mcsim`: Function `mcsim()` simulates metacommunity dynamics. By
    default, it produces a square-shaped landscape with randomly
    distributed habitat patches (x- and y-coordinates are drawn from a
    uniform distribution). If a distance matrix is given, the function
    simulates metacommunity dynamics in the specified landscape.
    Function `mcsim()` follows a general framework proposed by [Thompson
    et al. (2020)](https://doi.org/10.1111/ele.13568). However, it has
    several unique features that are detailed in [Terui et al
    (2021)](https://doi.org/10.1073/pnas.2105574118).

Mathematical details of these functions are described in [Terui et al
(2021)](https://doi.org/10.1073/pnas.2105574118).

# Citation

Current best citation for this package is:

Terui A, Kim S, Dolph CL, Kadoya T, Miyazaki Y. (2021) Emergent dual
scaling of riverine biodiversity. *Proceedings of the National Academy
of Sciences* 118: e2105574118

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
