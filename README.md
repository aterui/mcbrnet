mcbrnet
================

# Overview

`mcbrnet` is a collection of functions to perform metacommunity
simulation in branching networks. Main functions include:

-   `brnet()` generates a random branching network with the specified
    number of patches and probability of branching. The function returns
    adjacency and distance matrices, hypothetical
    environmental/disturbance values at each patch, and the number of
    patches upstream (akin to the watershed area in river networks). The
    output may be used in function `mcsim()` to simulate metacommunity
    dynamics in a branching network.

-   `mcsim()` simulates metacommunity dynamics. By default, it produces
    a square-shaped landscape with randomly distributed habitat patches
    (x- and y-coordinates are drawn from a uniform distribution). If a
    distance matrix is given, the function simulates metacommunity
    dynamics in the specified landscape. Function `mcsim()` follows a
    general framework proposed by [Thompson et
    al. (2020)](https://doi.org/10.1111/ele.13568). However, it has
    several unique features that are detailed in [Terui et al
    (2021)](https://doi.org/10.1073/pnas.2105574118).

-   `igpsim()` simulates three-species meta-food web dynamics with
    intraguild predation. This function shares many features with
    `mcsim()`. Currently under development.

-   `ggbrnet()` is a wrapper of `ggraph` functions for easy
    visualization of a network produced by `brnet()` .

-   `ptsource()` is a function to simulate propagation of environmental
    pollutants in a river network. Produced values may be plugged into
    `q` argument in `mcsim()`.

-   `adjtodist()` converts an adjacency matrix to a distance matrix.

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

# Change-log

#### v.1.3.0 (06/15/22)

-   add a new function `ptsource()`

-   add new arguments to `mcsim()`

-   update `ggbrnet()` to be compatible with piping

#### v.1.2.3 (04/13/22)

-   fix a bug in `fun_disp_mat()`

-   add `ggbrnet()`

-   add disturbance arguments to `mcsim()` (`p_disturb` & `m_disturb`)

#### v.1.2.2 (03/24/22)

-   fix a bug in `fun_igp()`

#### v.1.2.1 (03/09/22)

-   implement internal functions to `mcsim()` and `brnet()`

-   remove argument weighted_distance_matrix from `mcsim()` and
    `brnet()`

-   add argument dispersal_matrix to `mcsim()`

#### v.1.2.0 (03/09/22)

-   add a major function `igpsim()`

-   simplified `brnet()` and `mcsim()` by introducing internal
    sub-functions

#### v.1.1.1 (12/07/21)

-   add a local noise parameter for disturbance values to `brnet()`
    (argument `sd_disturb_lon`)

#### v.1.1.0 (08/02/21)

-   add disturbance arguments to `brnet()` added function `adjtodist()`

#### v.1.0.0 (05/03/21)

-   initial release

# Funding

This material is based upon work supported by the National Science
Foundation through the Division of Environmental Biology (DEB 2015634).
