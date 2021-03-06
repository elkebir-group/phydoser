% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phydose-mainfunctions.R
\name{CalcU}
\alias{CalcU}
\title{Generates the clonal prevalence matrix U by mutlplying the frequency matrix by
the inverse of the tree matrix}
\usage{
CalcU(treeMat, f)
}
\arguments{
\item{treeMat}{a binary matrix representing a tree with clones as rows and mutations as columns}

\item{f}{a frequency matrix with samples as rows and mutations as columns}
}
\value{
a matrix of clonal prevelance with samples as rows and clones as columns
}
\description{
Generates the clonal prevalence matrix U by mutlplying the frequency matrix by
the inverse of the tree matrix
}
