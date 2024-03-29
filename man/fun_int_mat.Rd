% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{fun_int_mat}
\alias{fun_int_mat}
\title{Interaction matrix}
\usage{
fun_int_mat(n_species, alpha, min_alpha, max_alpha, interaction_type)
}
\arguments{
\item{n_species}{Number of species in a metacommunity.}

\item{alpha}{Species interaction strength. Enabled if \code{interaction_type = "constant"}. Default \code{0}.}

\item{min_alpha}{Minimum value of a uniform distribution that generates species interaction strength. Enabled if \code{interaction_type = "random"}. Default \code{NULL}.}

\item{max_alpha}{Maximum value of a uniform distribution that generates species interaction strength. Enabled if \code{interaction_type = "random"}. Default \code{NULL}.}

\item{interaction_type}{\code{"constant"} or \code{"random"}. \code{"constant"} assumes the unique interaction strength of \code{alpha} for all pairs of species. \code{"random"} draws random numbers from a uniform distribution with \code{min_alpha} and \code{max_alpha}.}
}
\description{
Interaction matrix
}
