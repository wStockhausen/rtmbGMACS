% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calcGrowth.R
\name{calGrowth}
\alias{calGrowth}
\title{Calculate growth for all model categories across time}
\usage{
calGrowth(dims, info, params, verbose = FALSE, loopIC_ = TRUE)
}
\arguments{
\item{dims}{\itemize{
\item dimensions list
}}

\item{info}{\itemize{
\item info list (output list from \code{\link[=extractParamInfo_Growth]{extractParamInfo_Growth()}})
}}

\item{params}{\itemize{
\item RTMB parameters list with growth-specific elements
}}

\item{verbose}{\itemize{
\item flag to print diagnostic info
}}
}
\value{
TODO: might want to return a list of a list of matrices
}
\description{
Function to calculate growth for all model categories across time.
}
\details{
Growth matrix for any given y_, s_ is technically an upper triangle
block-diagonal, with non-zero elements only for z_row <= z_column (i.e., post-molt size)
AND {r,x,m,p}\emph{row == {r,x,m,p}}, followed by ??

At start, {r,x,m,p,z} has probability of molting prM(r,x,m,p,z), which splits
n_{r,x,m,p,z} into molting (mn_{r,x,m,p,z}) and non-molting (nn_{r,x,m,p,z}) components.
If terminal molt depends on pre-molt size, it should be evaluated now on molting animals
(e.g. immature->mature for mn).

The non-molting component should have p->max(p+1,p_max).

The molting component undergoes growth as a block-diagonal with non-zero transitions possible only
for z_row <= z_column (i.e., post-molt size) AND {r,x,m,p}_row == {r,x,m,p}_column, followed
by p->post-molt age 0.

If terminal molt depends on post-molt size, it would be evaluated now on molted crab
(e.g., immature-> mature).
}
