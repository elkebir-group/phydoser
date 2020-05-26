#' phydose example dataset of two trees
#' 
#' @format a list of four items:
#' \describe{
#'  \item{trees}{a list of a candidate set binary tree matrices with clones as rows and mutations as columns}
#'  \item{fmatrices}{a list of frequency matrices with samples as the rows and mutation frequencies as the columns}
#'  \item{u_list}{a list of clonal prevalence matrices with samples as the rows and clones as the columns}
#'  \item{graph}{a list of dot format string for each respective tree in the candate set}
#' }
#' 
"phydata"


#' phydose example distinguishing feature families 
#' @format a list of distinguishing feature families each containing a list of distinguishing f
#' eatures comprised of vectors of featurettes
#' 
"phydataDFF" 

#' phydose tutorial dataset of a patient (AML-38) with acute myeloid leukemia
#' 
#' @format a list of four items:
#' \describe{
#'  \item{trees}{a list of 316 candidate set binary tree matrices with clones as rows and mutations as columns}
#'  \item{fmatrices}{a  frequency matrix with samples as the rows and mutation frequencies as the columns}
#'  \item{u_list}{a list of 316 clonal prevalence matrices with samples as the rows and clones as the columns}
#'  \item{graph}{a list of 316 dot format strings for each respective tree in the candate set}
#' }
#' @source \url{https://www.biorxiv.org/content/10.1101/2020.02.07.925743v1.abstract}
"AML38"


#' AML38 corresponding distinguishing feature families 
#' @format a list of distinguishing feature families each containing a list of distinguishing f
#' eatures comprised of vectors of featurettes
#' 
"AML38DFF" 