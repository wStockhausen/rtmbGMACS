% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zzz_other_functions.R
\name{tibbleRowToString}
\alias{tibbleRowToString}
\title{Allow Quarto-based vignettes}
\usage{
tibbleRowToString(tbl_rw)
}
\arguments{
\item{tbl_rw}{- a row of a tibble}

\item{libname}{- library name (?)}

\item{pkgname}{- package name}
}
\value{
see [tools::vignetteEngine()].

a string
}
\description{
Function to allow Quarto-based vignettes

Function to convert a tibble row to a string.
}
\details{
Uses [tools::vignetteEngine()] to create a Quarto-based vignette engine
for this package.
}
