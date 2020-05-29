
<!-- README.md is generated from README.Rmd. Please edit that file -->

# phydoser

## Overview

phydoser is a package to aid in the design a single-cell sequencing
follow-up experiment from a previous bulk sequencing experiment.

  - `generateDistFeat` finds the distinguishing feature family of each
    tree in the input set  
  - `phydose` determines the number of cells to sequence in a follow-up
    single-cell experiment  
  - `CalcU` calculated the clonal prevalance matrix  
  - `plotTree` plots the tree along with any specified distinguishing
    features
  - `ReadJoint` reads in an input file into the data structures
    supported by phydose

## Installation

``` r
#development versions
devtools::install_github("elkebir-group/phydoser")
```

## Usage

``` r
library(phydoser)
```
