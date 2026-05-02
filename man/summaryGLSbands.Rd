\name{summaryGLSbands}
\alias{summaryGLSbands}
\title{Fit generalized least squares and plot confidence bands}
\description{
Estimate generalized least-squares regression coefficients and return a fitted plot with optional confidence bands.
}
\usage{
summaryGLSbands(OBJECT = NULL, Y = NULL, V_Y = NULL, PREDICTORS,
                INTERVALS = 20, CONTROLVARS = NULL, INTERCEPT = T,
                QUAD = F, BANDS = T, BAND_SIZE = 1, XLAB = "",
                YLAB = "", TITLELAB = "", XCOORDS = c(-2, 2),
                YCOORDS = c(0, 1))
}
\arguments{
  \item{OBJECT}{Optional object returned by \code{subSV}.}
  \item{Y}{Optional vector of summary parameters supplied directly.}
  \item{V_Y}{Optional sampling covariance matrix supplied directly.}
  \item{PREDICTORS}{Predictor vector or matrix for the GLS regression.}
  \item{INTERVALS}{Number of points used to evaluate confidence bands.}
  \item{CONTROLVARS}{Optional control variables appended to the design matrix.}
  \item{INTERCEPT}{Whether to include an intercept in the GLS model.}
  \item{QUAD}{Whether to include a quadratic predictor term.}
  \item{BANDS}{Whether to draw confidence bands.}
  \item{BAND_SIZE}{Multiplier applied to the confidence-band standard errors.}
  \item{XLAB}{X-axis label for the returned plot.}
  \item{YLAB}{Y-axis label for the returned plot.}
  \item{TITLELAB}{Plot title.}
  \item{XCOORDS}{X-axis limits.}
  \item{YCOORDS}{Y-axis limits.}
}
\value{
A plot object; the fitted coefficient table is also printed.
}
