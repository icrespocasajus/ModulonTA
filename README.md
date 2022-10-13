
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ModulonTA

## Title

Modulon Target Analysis

## Introduction

<!-- badges: start -->

[![R-CMD-check](https://github.com/icrespocasajus/ModulonTA/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/icrespocasajus/ModulonTA/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of ModulonTA is to analyze the downstream direct targets of
modulon constituent transcription factors in order to assess the regulon
similarity, redundancy and overlap between pairs of transcription
factors, or between transcription factors and a modulon regulatory core.
This information can be further used to identify modulon regulatory core
satellite transcription factors, or transcription factors that reinforce
the downstream regulatory effect of a given modulon regulatory core.

<img src="./man/figures/ModulonTA_flowchart.jpeg" width="100%" height="100%" style="display: block; margin: auto;" />

## Installation

You can install the development version of ModulonTA from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("icrespocasajus/ModulonTA")
```

## Example

This is a basic example which shows you how to perform a target analysis
on a specific modulon and plot the results:

``` r
library(ModulonTA)
library(corrplot)
#> corrplot 0.92 loaded

## Example
# For this example we will use the TILs dataset included in ModulonTA

network = network.TILs
modulons = modulons.TILs
modulon.query = '3'

results.target.analysis.modulon=target.analysis.modulon(net=network,mod=modulons,mod.query = modulon.query)
target.analysis.modulon.plot(data=results.target.analysis.modulon,feature = 'Redundancy')
```

<img src="man/figures/README-example-1.png" width="100%" style="display: block; margin: auto;" />

``` r

cc = cc.TILs
modulon.query = '3'
connected.component.query = 'cc.3'

results.target.analysis.modulon.wrt.cc.w.core = target.analysis.modulon.wrt.cc.manual.query.2(net = network,mod = modulons,cc = cc,mod.query = modulon.query,cc.query=connected.component.query)
```

``` r
satellites = Find.Sat(data = results.target.analysis.modulon.wrt.cc.w.core,feature = 'Redundancy',threshold = 0)
```

``` r
Discriminant.Analysis.data = DA.TILs
satellites.filtered = Filter.Sat(sat.data=satellites,DA.data = Discriminant.Analysis.data,DA=c("Any"),top.percent = 10)
```

## Author

Isaac Crespo

<img src="./man/figures/README-isaaccrespo_WEB_big.jpg" width="30%" />

Isaac Crespo, phD  
Senior Computational Scientist  
CHUV \| Department of Oncology \| George Coukos group  
Ludwig Institute for Cancer Research \| Lausanne Branch  
AGORA, Bugnon 25A, 1005 Lausanne, 4th floor, Room 026  
<isaaccrespo@hotmail.com>
