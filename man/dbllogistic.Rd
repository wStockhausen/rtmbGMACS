% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SelectivityFunctions.R
\name{dbllogistic}
\alias{dbllogistic}
\title{Calculate a logistic selectivity curve}
\usage{
dbllogistic(z, params, refZ = 0, debug = FALSE)
}
\arguments{
\item{z}{\itemize{
\item vector of sizes at which to compute selectivity curve
}}

\item{params}{\itemize{
\item 4-element vector of selectivity function parameters
}}

\item{refZ}{\itemize{
\item reference size
}}

\item{debug}{\itemize{
\item flag (T/F) to print debugging messages
}}
}
\value{
named vector with selectivity values at the elements of z
}
\description{
Function to calculate a logistic selectivity curve.
}
\details{

}
