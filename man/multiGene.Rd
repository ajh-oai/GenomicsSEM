\name{multiGene}
\alias{multiGene}
\title{Build multi-gene covariance structures for T-SEM}
\description{
Expand LDSC covariance and sampling covariance matrices to include multiple genes and their pairwise LD structure.
}
\usage{
multiGene(covstruc, Genes, LD, GeneSE = "F", Genelist = "F")
}
\arguments{
  \item{covstruc}{Output from multivariable LDSC.}
  \item{Genes}{Gene-level summary statistics containing beta and standard-error columns.}
  \item{LD}{Matrix of pairwise LD correlations among the genes.}
  \item{GeneSE}{Optional gene-variance standard error; the default uses the package convention.}
  \item{Genelist}{Optional vector of gene identifiers to retain from \code{Genes}.}
}
\value{
An expanded covariance structure containing multiple genes and their sampling covariance matrix.
}
