#' Cellwise regularized robust sparse regression
#'
#' @keywords package
#' @author Peng Su (peng.su@sydney.edu.au), Garth Tarr (garth.tarr@gmail.com), Samuel Muller (samuel.muller@mq.edu.au) and Suojin Wang (sjwang@stat.tamu.edu)
#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @importFrom robustbase Qn
#' @importFrom robustHD rlars sparseLTS robStandardize
#' @importFrom glmnet cv.glmnet
#' @importFrom mvtnorm rmvt
#' @importFrom stats coef median qchisq qt quantile rnorm
#' @importFrom utils data
#' @useDynLib regcell, .registration=TRUE
#' @name regcell
#' @references Su P, Tarr G, Muller S and Wang S (2023). CR-Lasso: Robust cellwise regularized sparse regression. arXiv:2307.05234
## usethis namespace: start
## usethis namespace: end
NULL


#' Bone mineral density data
#'
#' A subset of data from the European Bioinformatics Institute Array-Express repository.
#' The BMD data consists of gene expression measurements of 54,675 probes of 84 Norwegian women.
#'
#' @name datascreen
#' @docType data
#' @keywords datasets
#' @format A data frame with 84 rows and 101 columns.
#' @details Given the large number of variables in the dataset, a pre-screening step was implemented to identify the subset of variables that are most correlated with the outcome of interest, the total hip T-score.
#' To accomplish this, we first log-transformed all the predictors and then utilized the robust correlation estimate based on Winsorization.
#' The screened data comprise measurements of $p = 100$ genes from $n = 84$ Norwegian women.
#'
#' @source https://www.ebi.ac.uk/biostudies/files/E-MEXP-1618/Normarrayexpressdata.txt.magetab
#'  https://www.ebi.ac.uk/biostudies/files/E-MEXP-1618/E-MEXP-1618.sdrf.txt
#' @examples
#' data(datascreen)
NULL

