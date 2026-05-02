\name{addGenes}
\alias{addGenes}
\title{Combine LDSC and gene-level summary statistics for T-SEM}
\description{
Expand LDSC covariance and sampling covariance matrices to include one gene at a time for transcriptome-wide SEM analyses.
}
\usage{
addGenes(covstruc, Genes, GC = "standard")
}
\arguments{
  \item{covstruc}{Output from multivariable LDSC.}
  \item{Genes}{Gene-level summary statistics containing gene identifiers, heritability estimates, beta columns, and standard-error columns.}
  \item{GC}{Genomic-control correction to apply to standard errors: \code{"standard"}, \code{"conserv"}, or \code{"none"}.}
}
\value{
A list with expanded sampling covariance matrices, expanded covariance matrices, and gene identifiers.
}
