# README


-   [rtmbGMACS](#rtmbgmacs)
    -   [Installation](#installation)
    -   [Introduction](#introduction)
    -   [NOAA Disclaimer](#noaa-disclaimer)

<!--DO NOT VIEW THIS FILE USING THE VISUAL EDITOR!! (seems to screw things up)-->
<!-- README.md is generated from README.qmd. Please edit README.qmd, then render README.md using `quarto render README.qmd` in a terminal window. -->
<!-- use `quarto render README.qmd` in the terminal window to build README.md prior to committing to keep README.md up-to-date-->
<!-- don't forget to commit and push the resulting figure files, so they display on GitHub and CRAN.-->

# rtmbGMACS

<!-- badges: start -->
<!-- badges: end -->

rtmbGMACS is an [R](https://www.r-project.org/) package that, when
completed, will provide a RTMB-based implementation of the Generalized
Model for Assessing Crustacean Stocks
([GMACS](https://github.com/GMACS-project)) modeling framework that was
developed using AD Model Builder ([ADMB](http://www.admb-project.org)).
The ADMB version of GMACS is currently used to in stock assessment
models for a number of crab species in the Bering Sea and Aleutian
Islands to determine reference points for fishery management to the
North Pacific Fishery Management Council
([NPFMC](https://www.npfmc.org)) (see, for example, the stock assessment
reports for [Bristol Bay red king
crab](https://meetings.npfmc.org/CommentReview/DownloadFile?p=b98b90b2-88ab-43c2-9487-c12cdb4e0a25.pdf&fileName=BBRKC%20SAFE%202022%20Final.pdf)
and [snow
crab](https://meetings.npfmc.org/CommentReview/DownloadFile?p=fca55335-ad34-4896-9b1e-4c09aa8342ce.pdf&fileName=EBS%20Snow%20SAFE%20FINAL.pdf)).
RTMB ([R Bindings for Template Model
Builder](https://github.com/kaskr/RTMB)) is an R package that provides a
native R interface to Template Model Builder
([TMB](https://kaskr.github.io/adcomp)), an R package that allows the
user to develop and run C++-based statistical models (similar to ADMB)
using automatic differentiation in the R framework.

## Installation

**UNDER CONSTRUCTION: DON’T TRY THIS YET!!** You can install the current
(at present, development) version of rtmbGMACS from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("wStockhausen/rtmbGMACS")
```

## Introduction

The purpose of rtmbGMACS is to provide a population modeling and stock
assessment framework for, in particular, federally-managed
commercially-fished crab stocks in the US Bering Sea and Aleutian
Islands, although the software should be more widely applicable to (at
least) other crustacean stocks. Here, the modeling framework is referred
to as “gmacs” while the R package is referred to “rtmbGMACS”. The *R*
package will provide functionality to:

-   set up data inputs to a gmacs model
-   define associated likelihood functions
-   select functional forms for various population processes
-   define estimable parameters and associated priors or penalties
-   run the model to estimate parameters and stock status
-   determine required management-related quantities under different
    harvest strategies
-   create plots and tables from a fitted model
-   perform standard diagnostic analyses (e.g., retrospective analyses,
    likelihood profiling, MCMC)
-   and compare results from multiple models.

------------------------------------------------------------------------

## NOAA Disclaimer

This repository is a scientific product and is not an official
communication of the National Oceanic and Atmospheric Administration, or
the United States Department of Commerce. All NOAA GitHub project code
is provided on an ‘as is’ basis and the user assumes responsibility for
its use. Any claims against the Department of Commerce or Department of
Commerce bureaus stemming from the use of this GitHub project will be
governed by all applicable Federal law. Any reference to specific
commercial products, processes, or services by service mark, trademark,
manufacturer, or otherwise, does not constitute or imply their
endorsement, recommendation or favoring by the Department of Commerce.
The Department of Commerce seal and logo, or the seal and logo of a DOC
bureau, shall not be used in any manner to imply endorsement of any
commercial product or activity by DOC or the United States Government.

Software code created by U.S. Government employees is not subject to
copyright in the United States (17 U.S.C. §105). The United
States/Department of Commerce reserve all rights to seek and obtain
copyright protection in countries other than the United States for
Software authored in its entirety by the Department of Commerce. To this
end, the Department of Commerce hereby grants to Recipient a
royalty-free, nonexclusive license to use, copy, and create derivative
works of the Software outside of the United States.

------------------------------------------------------------------------

<img src="https://raw.githubusercontent.com/nmfs-general-modeling-tools/nmfspalette/main/man/figures/noaa-fisheries-rgb-2line-horizontal-small.png" height="75" alt="NOAA Fisheries">

[U.S. Department of Commerce](https://www.commerce.gov/) | [National
Oceanographic and Atmospheric Administration](https://www.noaa.gov) |
[NOAA Fisheries](https://www.fisheries.noaa.gov/)

------------------------------------------------------------------------

<figure>
<a
href="https://github.com/wStockhausen/rtmbGMACS/actions/workflows/pages/pages-build-deployment"><img
src="https://github.com/wStockhausen/rtmbGMACS/actions/workflows/pages/pages-build-deployment/badge.svg?branch=gh-pages"
alt="pages-build-deployment" /></a>
<figcaption>pages-build-deployment</figcaption>
</figure>
