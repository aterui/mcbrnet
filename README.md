mcbrnet: an R package for simulating metacommunity dynamics in a
branching network
================
Akira Terui
September 11, 2020

  - [Overview](#overview)
  - [Installation](#installation)
  - [Instruction](#instruction)
      - [`brnet()`](#brnet)
  - [References](#references)

# Overview

The package `mcbrnet` is composed of two functions: `brnet()` and
`mcsim()`.

  - `brnet`: Function `brnet()` generates a random branching network
    with the specified number of patches and probability of branching.
    The function returns adjacency and distance matrices, hypothetical
    environmental values at each patch, and the number of patches
    upstream (akin to the watershed area in river networks). The output
    may be used in function `mcsim()` to simulate metacommunity dynamics
    in a branching network.

  - `mcsim`: Function `mcsim()` simulates metacommunity dynamics. By
    default, it produces a square-shaped landscape with randomly
    distributed habitat patches (x- and y-coordinates are drawn from a
    uniform distribution). If a distance matrix is given, the function
    simulates metacommunity dynamics in the specified landscape.
    Function `mcsim()` follows a general framework proposed by Thompson
    et al. (2020). However, it has several unique features that are
    detailed in the following sections.

# Installation

The `mcbrnet` package can be installed with the following script:

``` r
#install.packages("remotes")
remotes::install_github("aterui/mcbrnet")
library(mcbrnet)
```

# Instruction

## `brnet()`

### Basic usage

The key arguments are the number of habitat patches (`n_patch`) and
probability of branching (`p_branch`), which users must specify. Using
these parameters, the function generates a branching network through the
following steps:

1.  Draw the number of branches in the network. An individual branch is
    defined as a series of connected patches from one confluence (or
    outlet) to the next confluence upstream (or upstream terminal). The
    number of branches in a network BR is drawn from a binomial
    distribution as BR \~ Binomial(N, P<sub>br</sub>), where N is the
    number of patches and P<sub>br</sub> is the branching probability.

2.  Draw the number of patches in each branch. The number of patches in
    each branch N<sub>br</sub> is drawn from a geometric distribution as
    N<sub>br</sub> \~ Ge(P<sub>br</sub>) but conditional on
    ΣN<sub>br</sub> = N.

3.  Organize branches into a bifurcating branching network.

The function returns:

  - `adjacency_matrix`: adjacency matrix.
  - `distance_matrix`: distance matrix. Distance between patches is
    measured as the number of patches required to reach from the focal
    patch to the target patch along the network.
  - `patch_df`: data frame. This data frame (`dplyr::tibble`) contains
    unique patch ID `patch_id`, unique branch ID `branch_id`,
    environmental value `environment` (see below for details), and the
    number of upstream contributing patches `n_patch_upstream`
    (including the focal patch itself; akin to the watershed area in
    river networks).

The following script produce a branching network with `n_patch = 50` and
`p_branch = 0.5`. By default, `brnet()` visualizes the generated network
using functions in packages `igraph` (Csardi and Nepusz 2006) and
`plotfunctions` (van Rij 2020) (`plot = FALSE` to disable):

``` r
net <- brnet(n_patch = 50, p_branch = 0.5)
```

![](README_files/figure-gfm/brnet_instruction_1-1.png)<!-- -->

Patches are colored in accordance with randomly generated environmental
values, and the size of patches is proportional to the number of patches
upstream. To view matrices, type the following script:

``` r
# adjacency matrix
# showing 10 patches for example
net$adjacency_matrix[1:10, 1:10]
```

    ##         patch1 patch2 patch3 patch4 patch5 patch6 patch7 patch8 patch9 patch10
    ## patch1       0      0      0      0      0      0      0      0      0       0
    ## patch2       0      0      1      0      0      0      0      0      0       0
    ## patch3       0      1      0      0      0      0      0      0      0       0
    ## patch4       0      0      0      0      1      0      0      0      0       0
    ## patch5       0      0      0      1      0      0      0      0      0       0
    ## patch6       0      0      0      0      0      0      1      0      0       0
    ## patch7       0      0      0      0      0      1      0      1      0       0
    ## patch8       0      0      0      0      0      0      1      0      1       0
    ## patch9       0      0      0      0      0      0      0      1      0       1
    ## patch10      0      0      0      0      0      0      0      0      1       0

``` r
# distance matrix
# showing 10 patches for example
net$distance_matrix[1:10, 1:10]
```

    ##         patch1 patch2 patch3 patch4 patch5 patch6 patch7 patch8 patch9 patch10
    ## patch1       0     10     11      4      5      3      4      5      6       7
    ## patch2      10      0      1     10     11      7      6      5      4       5
    ## patch3      11      1      0     11     12      8      7      6      5       6
    ## patch4       4     10     11      0      1      3      4      5      6       7
    ## patch5       5     11     12      1      0      4      5      6      7       8
    ## patch6       3      7      8      3      4      0      1      2      3       4
    ## patch7       4      6      7      4      5      1      0      1      2       3
    ## patch8       5      5      6      5      6      2      1      0      1       2
    ## patch9       6      4      5      6      7      3      2      1      0       1
    ## patch10      7      5      6      7      8      4      3      2      1       0

The following script lets you view branch ID, environmental values, and
the number of upstream contributing patches for each patch:

``` r
net$patch_df
```

    ## # A tibble: 50 x 4
    ##    patch_id branch_id environment n_patch_upstream
    ##       <int>     <dbl>       <dbl>            <dbl>
    ##  1        1         1      -0.481               50
    ##  2        2        12      -1.13                 2
    ##  3        3        12      -0.974                1
    ##  4        4        10       0.303                6
    ##  5        5        10       0.370                5
    ##  6        6         3      -0.431               26
    ##  7        7         3      -0.345               25
    ##  8        8         3      -0.595               24
    ##  9        9         3      -0.594               23
    ## 10       10        11      -0.180                6
    ## # ... with 40 more rows

### Customization

There is some flexibility in how to simulate environmental values, which
are determined through an autoregressive process, as detailed below:

1.  Environmental values for upstream terminal patches (i.e., patches
    with no upstream patch) are drawn from a uniform distribution as
    z<sub>1</sub> \~ Uniform(min<sub>env</sub>, max<sub>env</sub>).

2.  Downstream environmental values are determined by an autoregressive
    process as z<sub>x</sub> = ρz<sub>x-1</sub> + ε<sub>x</sub>, where
    ε<sub>x</sub> \~ Normal(0, σ<sup>2</sup><sub>env</sub>). At
    bifurcation patches (or confluence), the environmental value takes a
    weighted mean of the two contributing patches given the size of
    these patches *s* (the number of upstream contributing patches):
    z<sub>x</sub> = ω(ρz<sub>1,x-1</sub> + ε<sub>1,x</sub>) + (1 -
    ω)(ρz<sub>2,x-1</sub> + ε<sub>2,x</sub>), where ω =
    s<sub>1</sub>/(s<sub>1</sub> + s<sub>2</sub>).

Users may change the values of `min_env` (default: min<sub>env</sub> =
-1), `max_env` (max<sub>env</sub> = 1), `rho` (ρ = 1), and `sd_env`
(σ<sub>env</sub> = 0.1). Increasing the range of `min_env` and
`max_env` leads to greater variation in environmental values at upstream
terminals. `rho` (0 ≤ ρ ≤ 1) determines the strength of longitudinal
autocorrelation (the greater the stronger autocorrelation). `sd_env`
(σ<sub>env</sub> \> 0) determines the strength of local environmental
noise. The following script produces a network with greater
environmental variation at upstream terminals (z<sub>1</sub> \~
Uniform(-3, 3)), weaker longitudinal autocorrelation (ρ = 0.5), and
stronger local noises (σ<sub>env</sub> = 0.5).

``` r
net <- brnet(n_patch = 50, p_branch = 0.5,
             min_env = -3, max_env = 3, rho = 0.5, sd_env = 0.5)
```

![](README_files/figure-gfm/brnet_instruction_2-1.png)<!-- -->

# References

  - Csardi G, Nepusz T: The igraph software package for complex network
    research, InterJournal, Complex Systems 1695. 2006.
    <http://igraph.org>
  - Jacolien van Rij (2020). plotfunctions: Various Functions to
    Facilitate Visualization of Data and Analysis. R package version
    1.4. <https://CRAN.R-project.org/package=plotfunctions>
  - Thompson, P.L., Guzman, L.M., De Meester, L., Horváth, Z., Ptacnik,
    R., Vanschoenwinkel, B., Viana, D.S. and Chase, J.M. (2020), A
    process‐based metacommunity framework linking local and regional
    scale community ecology. Ecol Lett, 23: 1314-1329.
    <doi:10.1111/ele.13568>
