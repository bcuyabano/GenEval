% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GenEval.R
\name{popBreed}
\alias{popBreed}
\title{Simulates genotypes for generations of individuals}
\usage{
popBreed(M, pop.info)
}
\arguments{
\item{M}{A phased SNP matrix for generation zero.}

\item{pop.info}{The population info, in the same format generated by \code{simPopInfo}.}
}
\value{
A matrix with the SNPs for all individuals.
}
\description{
This function outputs SNP-genotypes for generations of individuals, based on the given SNPs for generation zero (founder population) and on the population info in the same format generated by \code{simPopInfo}.
}
\examples{
M <- simGeno(10,rep(0.5,3),phased=TRUE)
pop.info <- simPopInfo(c(10,15),1)
popBreed(M,pop.info)
}
\seealso{
\code{\link{simGeno}} \code{\link{simPopInfo}}
}
