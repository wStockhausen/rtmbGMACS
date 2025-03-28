% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SelectivityFunctions.R
\name{ascnormal3}
\alias{ascnormal3}
\title{Calculates a 2-parameter ascending normal selectivity curve}
\usage{
ascnormal3(z, params, refZ = 0, debug = FALSE)
}
\arguments{
\item{z}{\itemize{
\item dvector of sizes at which to compute function values
}}

\item{params}{\itemize{
\item dvar_vector of function parameters
}}

\item{refZ}{\itemize{
\item 2-element vector of reference sizes (see details)
}}
}
\value{
\itemize{
\item named vector of selectivity values
}
}
\description{
Function to calculate a 2-parameter ascending normal selectivity curve.
}
\details{
Calculates 2-parameter normal function parameterized by
\itemize{
\item params\link{1}: delta from max possible size (refZ\link{1}) at which ascending limb could reach 1
\item params\link{2}: selectivity at size=refZ\link{2}
}:
\code{refZ} is a 2-element vector with elements
\itemize{
\item refZ\link{1} - max possible size at which the curve could reach 1
\item refZ\link{2} - reference size at which curve reaches the value of param\link{2}
}
}
