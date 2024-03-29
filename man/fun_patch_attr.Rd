% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{fun_patch_attr}
\alias{fun_patch_attr}
\title{Weighted-mean patch attributes}
\usage{
fun_patch_attr(
  x,
  n_branch,
  mean_source,
  sd_source,
  sd_lon,
  m_distance,
  rho = 1,
  v_wa,
  logit = FALSE
)
}
\arguments{
\item{x}{Adjacency matrix to be converted}

\item{n_branch}{Number of branches}

\item{mean_source}{Mean value of source attribute}

\item{sd_source}{SD of source attribute}

\item{sd_lon}{SD of longitudinal change}

\item{m_distance}{Distance matrix}

\item{rho}{Longitudinal autocorrelation}

\item{v_wa}{Vector of watershed area}

\item{logit}{Logit transformation of mean source attribute}
}
\description{
Weighted-mean patch attributes
}
