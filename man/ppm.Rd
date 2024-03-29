% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcbrnet.R
\name{ppm}
\alias{ppm}
\title{Preferential prey model}
\usage{
ppm(n_species, n_basal, l, theta, cannibal = FALSE, lower_tri = TRUE)
}
\arguments{
\item{n_species}{Integer. Number of species}

\item{n_basal}{Integer. Number of basal species}

\item{l}{Interger. Expected number of links in the upper triangle}

\item{theta}{Numeric. Scale parameter of an exponential distribution. Smaller values indicate greater trophic specialization.}

\item{cannibal}{Logical. If \code{TRUE}, cannibalism allowed}

\item{lower_tri}{Logical. If \code{TRUE}, lower triangle elements of the matrix will be returned}
}
\description{
Preferential prey model
}
\author{
Akira Terui, \email{hanabi0111@gmail.com}
}
