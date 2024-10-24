% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SelectivityFunctions.R
\name{calcSelectivity}
\alias{calcSelectivity}
\title{Calculate values for a selectivity curve}
\usage{
calcSelectivity(type, z, params, ref = 0, debug = FALSE)
}
\arguments{
\item{type}{- the type of selectivity function to calculate}

\item{z}{- vector of values at which to calculate the function}

\item{params}{- the selectivity function parameters, as a vector}
}
\value{
vector matching size of z, with names given by elements of z
}
\description{
Calculate values for a selectivity curve.
}
\details{
"type" may be one of
\itemize{
\item{"sel_const"}
\item{"asclogistic"}
\item{'asclogistic5095'}
\item{'asclogistic50D95'}
\item{'dbllogistic'}
\item{'dbllogistic5095'}
\item{"ascnormal"}
\item{"ascnormal2"}
\item{"ascnormal2a"}
}
}
