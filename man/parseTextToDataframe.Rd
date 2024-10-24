% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MiscFunctions_ReadParse.R
\name{parseTextToDataframe}
\alias{parseTextToDataframe}
\title{Parse a text vector to create a dataframe}
\usage{
parseTextToDataframe(lns, iln, nlns, col_names = TRUE, show_col_types = FALSE)
}
\arguments{
\item{lns}{- character vector of text lines to parse}

\item{iln}{- index to line at which to start parsing}

\item{nlns}{- number of lines to parse}

\item{col_names}{- flag (T/F) or column names (see [readr::read_delim()])}

\item{show_col_types}{- flag (T/F) to show column types (see [readr::read_delim()])}
}
\value{
a tibble (see [tibble::tibble])
}
\description{
Parse a text vector to create a dataframe
}
\details{
Uses [readr::read_delim()] to parse the vector of text lines after
removing all trailing comments and replacing all white space with a
single space.
}
