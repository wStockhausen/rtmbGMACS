% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SelectivityFunctions.R
\name{asclogistic5095}
\alias{asclogistic5095}
\title{Calculate an ascending logistic function parameterized by z50 and z95}
\usage{
asclogistic5095(z, params, refZ = 0, debug = FALSE)
}
\arguments{
\item{z}{- vector of sizes at which to compute selectivities}

\item{params}{- 2-element vector with selectivity function parameters}

\item{refZ}{- reference size}

\item{debug}{- flag (T/F) to print debugging messages}
}
\value{
named vector with selectivity values at the elements of z
}
\description{
Function to calculate an ascending logistic function parameterized by z50 and z95.
}
\details{
The parameter values are
\itemize{
\item{params[1] - z50: size at which selectivity = 0.50 (logit-scale mean)}
\item{params[2] - z95: size at which selectivity = 0.95}
}

If `refZ`>0, `refZ`=fully-selected size.
If `refZ`<0, function is normalized to max.
If `refZ`=0, no re-scaling is done.
}
