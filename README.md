
<!-- README.md is generated from README.Rmd. Please edit that file -->

# phydoser

## Overview

phydoser is a package to aid in the design a single-cell sequencing
follow-up experiment from a previous bulk sequencing experiment. The
main functions of phydoser include:

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

To install phydoser, you will first need to install the devtools
package:

``` r
install.packages("devtools")
```

Then you can easily install phydoser from Github:

``` r
#development version
devtools::install_github("elkebir-group/phydoser")
```

## Usage

``` r
library(phydoser)
```

phydoser comes with three test datasets that you can use to explore the
packages capabilities: `AML38`, `phydata`, and `multsample`.
