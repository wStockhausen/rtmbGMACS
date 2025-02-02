% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SelectivityFunctions.R
\name{dblnormal4}
\alias{dblnormal4}
\title{Calculates a 4-parameter normal selectivity curve}
\usage{
dblnormal4(z, params, refZ = 0, debug = FALSE)
}
\arguments{
\item{z}{- dvector of sizes at which to compute function values}

\item{params}{- dvar_vector of function parameters}

\item{refZ}{- size at which function = 1 (i.e., fully-selected size) [double] NOTE: ignored!}
}
\value{
- named vector of selectivity values
}
\description{
Function to calculate a 4-parameter normal selectivity curve.
}
\details{
Calculates 4-parameter normal function parameterized by
\itemize{
     \item params[1]: size at which ascending limb reaches 1
     \item params[2]: width of ascending limb
     \item params[3]: offset to size at which descending limb departs from 1
     \item params[4]: width of descending limb
}
}
