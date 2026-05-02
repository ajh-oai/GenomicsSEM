\name{summaryGLS}
\alias{summaryGLS}
\title{Fit generalized least squares to summary parameters}
\description{
Estimate generalized least-squares regression coefficients for a vector of summary parameters and their sampling covariance matrix.
}
\usage{
summaryGLS(OBJECT = NULL, Y = NULL, V_Y = NULL, PREDICTORS,
           INTERCEPT = T)
}
\arguments{
  \item{OBJECT}{Optional object returned by \code{subSV}.}
  \item{Y}{Optional vector of summary parameters supplied directly.}
  \item{V_Y}{Optional sampling covariance matrix supplied directly.}
  \item{PREDICTORS}{Predictor vector or matrix for the GLS regression.}
  \item{INTERCEPT}{Whether to include an intercept in the GLS model.}
}
\value{
Prints a coefficient table containing estimates, p-values, standard errors, and Z statistics.
}
