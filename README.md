mcbrnet: An R package for simulating metacommunity dynamics in a
branching network
================
Akira Terui
9/10/2020

  - [Overview](#overview)
  - [Installation](#installation)
  - [Instruction](#instruction)
      - [`brnet()`](#brnet)
      - [`brnet()`: Customization](#brnet-customization)
  - [References](#references)

# Overview

The package `mcbrnet` is composed of two functions: `brnet()` and
`mcsim()`.

  - `brnet`: Function `brnet()` generates a random branching network
    with the specified number of patches and probability of branching.
    The function returns adjacency and distance matrices, hypothetical
    environmental values at each patch, and the number of patches
    upstream (akin to watershed area in river networks). The output may
    be used in function `mcsim()` to simulate metacommunity dynamics in
    a branching network.

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

### Basic function

The function `brnet()` generates a random branching network. The key
arguments are the number of habitat patches (`n_patch`) and probability
of branching (`p_branch`), which the user must specify (otherwise, it
will be set to be `n_patch = 100` and `p_branch = 0.5`). The branching
network will be generated through the following steps:

1.  Determine the number of branches in the network. An individual
    branch is defined as a series of connected patches from one
    confluence (or outlet) to the next confluence upstream (or upstream
    terminal). The number of branches in a network `n_branch` is drawn
    from a binomial distribution as `n_branch = rbinom(n = 1, size =
    n_patch, prob = p_branch)`.

2.  Determine the number of patches in each branch. The number of
    patches in each branch `v_n_patch_branch` is drawn from a geometric
    distribution as `v_n_patch_branch = rgeom(n = n_branch, prob =
    p_branch) + 1`.

3.  Organize branches into a bifurcating branching network.

The following script produce a branching network with `n_patch = 50` and
`p_branch = 0.5`, returning adjacency and distance matrices. By default,
the function visualizes the generated network using functions in package
`igraph` (Csardi and Nepusz 2006):

``` r
net <- brnet(n_patch = 50, p_branch = 0.5)
```

![](README_files/figure-gfm/brnet_instruction_1-1.png)<!-- -->

To view matrices, type the following script:

``` r
# adjacency matrix
# select initial 10 patches for example
net$adjacency_matrix[1:10, 1:10]
```

    ##       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
    ##  [1,]    0    0    0    0    0    0    0    0    0     0
    ##  [2,]    0    0    1    0    0    0    0    0    0     0
    ##  [3,]    0    1    0    0    0    0    0    0    0     0
    ##  [4,]    0    0    0    0    1    0    0    0    0     0
    ##  [5,]    0    0    0    1    0    0    0    0    0     0
    ##  [6,]    0    0    0    0    0    0    1    0    0     0
    ##  [7,]    0    0    0    0    0    1    0    1    0     0
    ##  [8,]    0    0    0    0    0    0    1    0    1     0
    ##  [9,]    0    0    0    0    0    0    0    1    0     1
    ## [10,]    0    0    0    0    0    0    0    0    1     0

``` r
# distance matrix
# select initial 10 patches for example
net$distance_matrix[1:10, 1:10]
```

    ##       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
    ##  [1,]    0   10   11    4    5    3    4    5    6     7
    ##  [2,]   10    0    1   10   11    7    6    5    4     5
    ##  [3,]   11    1    0   11   12    8    7    6    5     6
    ##  [4,]    4   10   11    0    1    3    4    5    6     7
    ##  [5,]    5   11   12    1    0    4    5    6    7     8
    ##  [6,]    3    7    8    3    4    0    1    2    3     4
    ##  [7,]    4    6    7    4    5    1    0    1    2     3
    ##  [8,]    5    5    6    5    6    2    1    0    1     2
    ##  [9,]    6    4    5    6    7    3    2    1    0     1
    ## [10,]    7    5    6    7    8    4    3    2    1     0

## `brnet()`: Customization

# References

  - Csardi G, Nepusz T: The igraph software package for complex network
    research, InterJournal, Complex Systems 1695. 2006.
    <http://igraph.org>
  - Thompson, P.L., Guzman, L.M., De Meester, L., Horváth, Z., Ptacnik,
    R., Vanschoenwinkel, B., Viana, D.S. and Chase, J.M. (2020), A
    process‐based metacommunity framework linking local and regional
    scale community ecology. Ecol Lett, 23: 1314-1329.
    <doi:10.1111/ele.13568>
