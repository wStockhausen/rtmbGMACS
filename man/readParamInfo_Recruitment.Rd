% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readParamInfo_Recruitment.R
\name{readParamInfo_Recruitment}
\alias{readParamInfo_Recruitment}
\title{Read recruitment-related parameter info}
\usage{
readParamInfo_Recruitment(conn, debug = FALSE)
}
\arguments{
\item{conn}{\itemize{
\item connection (or filename) to read from
}}

\item{debug}{\itemize{
\item flag to print extra info
}}
}
\value{
nested list with elements (see \strong{Details})
}
\description{
Function to read recruitment-related parameter info.
}
\details{
The returned list has elements
\itemize{
\item{lstPrcs - process specification list, has elements}
\itemize{
\item{n - number of specifications}
\item{dfr - dataframe with process info}
}
\item{lstTBs - process specification list, has elements}
\itemize{
\item{n - number of time blocks}
\item{dfr - dataframe with time blocks info}
}
\item{lstTotRec - total recruitmnet specification list, has elements}
\item{lstXR - recruitment sex ratio list, has elements}
\item{lstRecZD - recruitment size distribution list, has elements}
}
where \code{lstTotRec}, \code{lstXR}, and \code{lstRecZD} are each lists having elements
\itemize{
\item{nFcns - number of functions specified}
\item{dfrFcns - dataframe with function info}
\item{nPars - number of main parameters specified}
\item{dfrPars - dataframe with main parameters info}
\item{nRWs - number of random walk parameters specified}
\item{dfrRWs - dataframe with random walk parameters info}
\item{nECs - number of environmental covariates specified}
\item{dfrECs - dataframe with environmental covariates info}
}
}
