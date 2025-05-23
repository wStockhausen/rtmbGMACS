% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calcWatZ.R
\name{calcWatZ_data}
\alias{calcWatZ_data}
\title{"Calculate" weight-at-size based on "data" inputs}
\usage{
calcWatZ_data(dims, dfrDims2Pars, par_vals)
}
\arguments{
\item{dims}{\itemize{
\item model dimensions tibble at which to "calculate" weight-at-size
}}

\item{dfrDims2Pars}{\itemize{
\item tibble map from model dims to parameter indices
}}

\item{par_vals}{\itemize{
\item vector of potential parameter (i.e., weight-at-size) values
}}
}
\value{
vector of weights-at-size for the given model dimensions in \code{dims}
}
\description{
Function to "calculate" weight-at-size based on "data" inputs.
}
\details{
To use this function, values for weight-at-size (i.e., "data") have
been given and expanded to all relevant model dimensions. Consequently, the function
merely looks up (hence the quotes around "calculate" in the title and description)
and returns the values associated with the given \code{dims} tibble.

Sample code

\if{html}{\out{<div class="sourceCode">}}\preformatted{dirPrj = rstudioapi::getActiveProject();
conn=file.path(dirPrj,"testing/testAllometry/inputSpecs_Allometry.data.txt");
res = readParamInfo_Allometry(conn,verbose=FALSE);
res1 = extractParameters_Allometry(res,dims$dms_yrxmaz,verbose=FALSE);
sel_dims = (dims$dms_yrxmaz |> dplyr::filter(y==2021))[1:3,];#--just extract values at 3 dimension combinations
sel_pars = calcWatZ_data(sel_dims,res1$dfrDims2Pars,res1$params);
dplyr::bind_cols(sel_dims,par_val=sel_pars)
# A tibble: 3 × 8
#     sparse_idx y    r     x      m        a         z       par_vals
#        <int> <fct> <fct> <fct>  <fct>    <fct>     <fct>      <dbl>
# 1      11777 2021  EBS   female immature new_shell 27    0.00000605
# 2      11778 2021  EBS   female immature new_shell 32    0.00000976
# 3      11779 2021  EBS   female immature new_shell 37    0.0000147
}\if{html}{\out{</div>}}
}
