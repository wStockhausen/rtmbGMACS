% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MiscFunctions_ReadParse.R
\name{parseLines}
\alias{parseLines}
\title{Parse multiple text lines into a list with named elements}
\usage{
parseLines(lns, names)
}
\arguments{
\item{lns}{\itemize{
\item the text lines parse
}}

\item{names}{\itemize{
\item the names of the elements of the output list
}}
}
\value{
named list of parsed values
}
\description{
Function to parse multiple text lines into a list with named elements.
}
\details{
The returned list has element names matching the \code{names} provided.
Each non-white space line that doesn't start with a comment character ('#') will
result in an element of the returned list. Vector-valued elements should have
all values on the same line. Matrices and other two-dimensional (or higher)
constructs cannot be paresd corretly using this function.

IMPORTANT: A counter with the name \code{iln} must be defined in a parent environment
to this function. The parsing of lines in \code{lns} uses the value of iln
at the time the function is called as the starting index in \code{lns}.
When the function exits, the value of \code{iln} indicates the index of the last element
examined.
}
