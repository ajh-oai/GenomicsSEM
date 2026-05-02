\name{subSV}
\alias{subSV}
\title{Subset LDSC covariance and sampling covariance elements}
\description{
Extract selected covariance or correlation elements and their sampling covariance submatrix from LDSC output.
}
\usage{
subSV(LDSC_OBJECT = NULL, SMATRIX = NULL, VMATRIX = NULL, INDEXVALS,
      TYPE = "S")
}
\arguments{
  \item{LDSC_OBJECT}{Optional LDSC output object.}
  \item{SMATRIX}{Optional covariance or correlation matrix supplied directly.}
  \item{VMATRIX}{Optional sampling covariance matrix supplied directly.}
  \item{INDEXVALS}{Indices of the lower-triangular elements to retain.}
  \item{TYPE}{Matrix type to extract from \code{LDSC_OBJECT}: \code{"S"}, \code{"S_Stand"}, or \code{"R"}.}
}
\value{
A list containing the subsetted sampling covariance matrix and covariance vector.
}
