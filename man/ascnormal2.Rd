% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SelectivityFunctions.R
\name{ascnormal2}
\alias{ascnormal2}
\title{Calculate ascending normal selectivity curve}
\usage{
ascnormal2(z, params, refZ = 0, debug = FALSE)
}
\arguments{
\item{z}{\itemize{
\item vector of sizes at which to compute selectivity values
}}

\item{params}{\itemize{
\item 2-element vector of selectivity function parameters
}}

\item{refZ}{\itemize{
\item size at which function reaches params\link{2}
}}

\item{debug}{\itemize{
\item flag (T/F) to print debugging messages
}}
}
\value{
named vector with selectivity values at the elements of z
}
\description{
Function to calculate ascending normal selectivity curve.
}
\details{
The parameter vector has values
\itemize{
\item{params\link{1}: selectivity at size = refZ}
\item{params\link{2}: size at which ascending limb reaches 1}
}
}
