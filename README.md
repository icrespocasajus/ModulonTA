
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ModulonTA

## Title

Modulon Target Analysis

## Introduction

The goal of ModulonTA is to analyze the downstream direct targets of
modulon constituent transcription factors in order to assess the regulon
similarity, redundancy and overlap between pairs of transcription
factors, or between transcription factors and a modulon regulatory core.
This information can be further used to identify modulon regulatory core
satellite transcription factors, or transcription factors that reinforce
the downstream regulatory effect of a given modulon regulatory core.

<img src="./man/figures/ModulonTA_flowchart.jpg" width="100%" height="100%" style="display: block; margin: auto;" />

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
feature = 'Redundancy'

results.target.analysis.modulon=target.analysis.modulon(net=network,mod=modulons,mod.query = modulon.query)
target.analysis.modulon.plot(data=results.target.analysis.modulon,feature = feature)
```

<img src="man/figures/README-example-1.png" width="100%" style="display: block; margin: auto;" />

In order to identify the modulon regulatory core satellites we are going
to perform the modulon target analysis with respect the query modulon
regulatory core as a whole:

``` r


network = network.TILs
modulons = modulons.TILs
connected.components = cc.TILs
modulon.query = '3'
connected.component.query = 'cc.3'

results.target.analysis.modulon.wrt.cc.w.core = target.analysis.modulon.wrt.cc.manual.query(net = network,mod = modulons,cc = connected.components,mod.query = modulon.query,cc.query=connected.component.query)
head(results.target.analysis.modulon.wrt.cc.w.core[[1]][['Redundancy']])
#>            TF Redundancy
#> Crem     Crem 0.31269350
#> Fli1     Fli1 0.19814241
#> Tbx21   Tbx21 0.16099071
#> Runx3   Runx3 0.12074303
#> Gabpb1 Gabpb1 0.11145511
#> Maf       Maf 0.05882353
```

The object results.target.analysis.modulon.wrt.cc.w.core contains target
analysis results using 3 different metrics: ‘Similarity’, ‘Redundancy’
and ‘Overlap’.

Using the results of the target analysis with respect to the modulon
regulatory core we are going to search for the satellites:

``` r
regulatory.core.elements = connected.components[[modulon.query]][[connected.component.query]]
satellites = Find.Sat(data = results.target.analysis.modulon.wrt.cc.w.core,feature = "Redundancy",threshold = 0,core = regulatory.core.elements )
satellites
#> $`3__cc.3`
#>  [1] "Zmiz1"   "Smarca4" "Eomes"   "Runx2"   "Rora"    "Nr3c1"   "Cebpd"  
#>  [8] "Gata3"   "Tead1"   "Batf"    "Foxn2"   "Srebf2"  "Zfp821"  "Srf"    
#> [15] "Stat5a"  "Tbp"     "Mxd3"    "Vezf1"   "Psmd12"
```

We can filter out those satellites with low discriminant power and
retain only those within top 10% of the discriminant power ranking for
any of the T cell subpopulations:

``` r
Discriminant.Analysis.data = DA.TILs
satellites.filtered = Filter.Sat(sat.data=satellites,DA.data = Discriminant.Analysis.data,DA=c("Any"),top.percent = 10)
satellites.filtered
#> $`3__cc.3`
#>  [1] "Zmiz1"   "Smarca4" "Eomes"   "Runx2"   "Rora"    "Nr3c1"   "Cebpd"  
#>  [8] "Gata3"   "Tead1"   "Batf"    "Foxn2"   "Srebf2"  "Zfp821"  "Srf"    
#> [15] "Mxd3"
```

We can plot a summmary heatmap including the modulon regulatory core and
satellites showing the pair-wise target redundancy, redundancy with
respect the core and the discriminant power across T cell subpopulations
from the OPLS-DA:

``` r
RegAUC = RegAUC.TILs
target.analysis.heatmap(
            net=network,
            mod=modulons,
            cc=connected.components,
            mod.query =  modulon.query,
            cc.query = connected.component.query,
            feature='Redundancy',
            sat=satellites.filtered,
            DA.data=Discriminant.Analysis.data,
            DA='Any',
            RegAUC=RegAUC,
            results.target.analysis.modulon=results.target.analysis.modulon,
            color= 'YlGn')
```

<img src="man/figures/README-heatmap-1.png" width="100%" height="100%" />

The annotation barplot on the right hand corresponds to the target
analysis with respect to the modulon regulatory core as a whole; this
value is the one used to select the satellites. The Core.Membership
annotation corresponds with the correlation of the regulon activity and
the firs principal component (PC1) of a PCA including only the modulon
regulatory core constituent elements; high membership might point out
actual modulon regulatory core elements wrongly excluded from the the
core.

## Author

Isaac Crespo

<img src="./man/figures/README-isaaccrespo_WEB_big.jpg" width="30%" />

Isaac Crespo, phD  
Senior Computational Scientist  
CHUV \| Department of Oncology \| George Coukos group  
Ludwig Institute for Cancer Research \| Lausanne Branch  
AGORA, Bugnon 25A, 1005 Lausanne, 4th floor, Room 026  
<isaaccrespo@hotmail.com>
