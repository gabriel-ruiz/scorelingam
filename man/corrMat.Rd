% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrapper-functions.r
\name{corrMat}
\alias{corrMat}
\title{Quickly estimate a Pearson correlation matrix using RcppArmadillo.}
\usage{
corrMat(X)
}
\arguments{
\item{X}{An n by p data matrix which is used to estimate the correlation matrix.}
}
\value{
A p by p estimate of the correlation matrix.
}
\description{
Quickly estimate a Pearson correlation matrix using RcppArmadillo.
}
