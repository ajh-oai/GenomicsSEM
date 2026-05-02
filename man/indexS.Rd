\name{indexS}
\alias{indexS}
\title{Index unique elements of a covariance or correlation matrix}
\description{
Create a symmetric index matrix for the unique lower-triangular elements of an LDSC covariance or correlation matrix.
}
\usage{
indexS(LDSC_OBJECT = NULL, MATRIX = NULL, R = F)
}
\arguments{
  \item{LDSC_OBJECT}{Optional LDSC output object containing an \code{S} matrix.}
  \item{MATRIX}{Optional covariance or correlation matrix supplied directly.}
  \item{R}{Whether to index only off-diagonal correlation elements.}
}
\value{
A symmetric integer matrix mapping matrix cells to lower-triangular element indices.
}
