#'
#' @title Allow Quarto-based vignettes
#'
#' @description Function to allow Quarto-based vignettes
#'
#' @param libname - library name (?)
#' @param pkgname - package name
#'
#' @return see [tools::vignetteEngine()].
#'
#' @details Uses [tools::vignetteEngine()] to create a Quarto-based vignette engine
#' for this package.
#'
#' @importFrom tools vignetteEngine
#'
# .onLoad <- function(libname, pkgname) {
#   tools::vignetteEngine(name = "quarto",
#                         package = "rtmbGMACS",
#                         pattern = "[.]qmd$",
#                         weave = function(file, ..., encoding = "UTF-8") {
#                           message("Hello from custom Quarto vignette engine!")
#                           ## TODO: What is Quarto's expectation about encoding?
#                           ## NB: output_format = "all" below might make a better default
#                           quarto::quarto_render(file, ..., output_format = "all")
#                         },
#                         tangle = tools::vignetteEngine("knitr::rmarkdown")$tangle,
#                         aspell = tools::vignetteEngine("knitr::rmarkdown")$aspell
#                         )
# }

#'
#' @title Convert a tibble row to a string
#'
#' @description Function to convert a tibble row to a string.
#'
#' @param tbl_rw - a row of a tibble
#'
#' @return a string
#'
#' @export
#'
tibbleRowToString<-function(tbl_rw){
  str = "";
  for (cn in rev(names(tbl_rw))){
    if (is.factor(tbl_rw[[cn]])){
      str = paste0(as.character(tbl_rw[[cn]])," ",str);
    } else {
      str = paste0(as.character(tbl_rw[[cn]])," ",str);
    }
  }
  return(str)
}

#--tibbleRowToString(tbl_rw)
