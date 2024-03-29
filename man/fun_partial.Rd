% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{fun_partial}
\alias{fun_partial}
\title{Yield partial derivatives}
\usage{
fun_partial(r, a, i, x0, model)
}
\arguments{
\item{r}{Numeric. Intrinsic growth rate of species i.}

\item{a}{Numeric vector. Interaction coefficients.}

\item{i}{Integer. Species numeric ID.}

\item{x0}{Numeric vector. Specify a vector of equilibrium densities for modeled species.}

\item{model}{Character. Specify a model type. Either \code{"ricker"} (Ricker), \code{"bh"} (Beverton-Holt), or \code{"glv"} (Generalized Lotka-Volterra).}
}
\description{
Yield partial derivatives
}
