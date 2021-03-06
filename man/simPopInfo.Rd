% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GenEval.R
\name{simPopInfo}
\alias{simPopInfo}
\title{Simulate a population info}
\usage{
simPopInfo(
  n,
  n.gen = 0,
  p.male = 0.5,
  pop.info = NULL,
  pedigree = FALSE,
  pre.id = ""
)
}
\arguments{
\item{n}{The size of each generation simulated. If a single number is provided, all generations will be of the same size.}

\item{n.gen}{The number of generations to simulate after generation zero. The generation counter starts at zero, being that the founder population. By default, \code{n.gen=0}, and if \code{n.gen} is set to any positive number, e.g. \code{n.gen=1}, that means that the info on the founder generation and on one generation after that will be generated.}

\item{p.male}{The proportion of male individuals in the population.}

\item{pop.info}{The information on a base population on which the code should build new generations.}

\item{pedigree}{Should the pedigree be generated? Default is \code{pedigree=FALSE}.}

\item{pre.id}{An ID label to paste at the begining of the default ID generated by the function.}
}
\value{
A data frame with the population info, and the pedigree when required.
\describe{
\item{info}{The population info.}
\item{ped}{The pedigree matrix.}
}
}
\description{
This function simulates a population info, with individuals' ID, sex, generation and parent information.
}
\examples{
simPopInfo(10,0)
simPopInfo(10,2)
simPopInfo(3,1,pedigree=TRUE)
simPopInfo(10,2,0.25)
simPopInfo(c(5,10,15),2)
}
