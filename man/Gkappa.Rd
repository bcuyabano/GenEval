% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GenEval.R
\name{Gkappa}
\alias{Gkappa}
\title{Kappa statistics}
\usage{
Gkappa(G1, G2 = NULL, U2 = NULL, L2 = NULL, Kmat = FALSE)
}
\arguments{
\item{G1}{The Reference matrix.}

\item{G2}{The matrix which we want to compare with G1.}

\item{U2}{The eigen-vectors of G1. If not provided, the function will obtain them internally.}

\item{L2}{The eigen-values of G2. If not provided, the function will obtain them internally.}

\item{Kmat}{Should the entire Kappa matrix be returned? Default is \code{Kmat=FALSE}.}
}
\value{
\describe{
\item{lambda}{The eigen-values of G2.}
\item{kappa}{The kappa statistics that compare G2 to G1.}
\item{Kmat}{The whole Kappa matrix.}
}
}
\description{
This function calculates the kappa statistics, to compare two variance-covariance matrices (G1 and G2) defined for the same group of individuals.
Kappa is computed using the eigen-decomposition of the matrices, \code{G=ULU'} in which U and L are the eigen-vectors and eigen-values, respectively. \code{Kmat = U2'(G1)U2}, and \code{kappa=diag(Kmat)}.
Finally, to compare G1 and G2, it is enough to plot \code{lambda} against \code{kappa}. The more linear their relationship, the more equivalent matrices G1 and G2 are, in terms of variance structure.
}
\examples{
# simulate uncorrelated genotypes
M <- simGeno(n=100,AlleleFreq=rep(0.5,200))
# three GRMs from the simulated genotypes
# the first with 100 SNPs
G1 <- mkGRM(M[,1:100])
# the second with 100 SNPs that don't overlap those from the first
G2 <- mkGRM(M[,101:200])
# the third with all SNPs (100 SNPs overlap with those from the first matrix)
G3 <- mkGRM(M)
#-------------------------------------------------------------------#
# compare GRMs using kappa
# verify how that is different from directly comparing the values
# verify how that is also different from comparing the eigen-values
#-------------------------------------------------------------------#
###########
# G1 x G2 #
###########
par(mar=c(4.5,4.5,1,1))
# plot G1 x G2
plot(c(diag(G1),G1[upper.tri(G1)]),c(diag(G2),G2[upper.tri(G2)]),xlab="G1",ylab="G2")
abline(0,1,lty=2,col=2)
# plot eigen-values(G1) x eigen-values(G2)
plot(eigen(G1)$values,eigen(G2)$values,xlab="eigen-values(G1)",ylab="eigen-values(G2)")
abline(0,1,lty=2,col=2)
# plot eigen-values(G2)xkappa
tmp <- Gkappa(G1,G2)
plot(tmp$lambda,tmp$kappa,xlim=range(tmp),ylim=range(tmp),xlab="lambda",ylab="kappa")
abline(0,1,lty=2,col=2)
###########
# G1 x G3 #
###########
par(mar=c(4.5,4.5,1,1))
# plot G1xG3
plot(c(diag(G1),G1[upper.tri(G1)]),c(diag(G3),G3[upper.tri(G3)]),xlab="G1",ylab="G3")
abline(0,1,lty=2,col=2)
# plot eigen-values(G1)xeigen-values(G3)
plot(eigen(G1)$values,eigen(G3)$values,xlab="eigen-values(G1)",ylab="eigen-values(G3)")
abline(0,1,lty=2,col=2)
# plot eigen-values(G3)xkappa
tmp <- Gkappa(G1,G3)
plot(tmp$lambda,tmp$kappa,xlim=range(tmp),ylim=range(tmp),xlab="lambda",ylab="kappa")
abline(0,1,lty=2,col=2)

# simulate correlated genotypes
M <- simGeno(n=100,AlleleFreq=rep(0.5,1000),LD=TRUE,phased=FALSE)
# three GRMs from the simulated genotypes
# the first with 500 SNPs
G1 <- mkGRM(M[,1:500])
# the second with 500 SNPs that don't overlap those from the first
G2 <- mkGRM(M[,501:1000])
# the third with all SNPs (500 SNPs overlap with those from the first matrix)
G3 <- mkGRM(M)
#-------------------------------------------------------------------#
# compare GRMs using kappa
# verify how that is different from directly comparing the values
# verify how that is also different from comparing the eigen-values
#-------------------------------------------------------------------#
###########
# G1 x G2 #
###########
par(mfrow=c(2,2),mar=c(4.5,4.5,1,1))
# plot G1xG2
plot(c(diag(G1),G1[upper.tri(G1)]),c(diag(G2),G2[upper.tri(G2)]),xlab="G1",ylab="G2")
abline(0,1,lty=2,col=2)
# plot eigen-values(G1)xeigen-values(G2)
plot(eigen(G1)$values,eigen(G2)$values,xlab="lambda(G1)",ylab="lambda(G2)")
abline(0,1,lty=2,col=2)
# plot eigen-values(G2)xkappa
tmp <- Gkappa(G1,G2)
plot(tmp$lambda,tmp$kappa,xlim=range(tmp),ylim=range(tmp),xlab="lambda",ylab="kappa")
abline(0,1,lty=2,col=2)
plot(NA,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",axes=FALSE)
###########
# G1 x G3 #
###########
par(mfrow=c(2,2),mar=c(4.5,4.5,1,1))
# plot G1xG3
plot(c(diag(G1),G1[upper.tri(G1)]),c(diag(G3),G3[upper.tri(G3)]),xlab="G1",ylab="G3")
abline(0,1,lty=2,col=2)
# plot eigen-values(G1)xeigen-values(G3)
plot(eigen(G1)$values,eigen(G3)$values,xlab="lambda(G1)",ylab="lambda(G3)")
abline(0,1,lty=2,col=2)
# plot eigen-values(G3)xkappa
tmp <- Gkappa(G1,G3)
plot(tmp$lambda,tmp$kappa,xlim=range(tmp),ylim=range(tmp),xlab="lambda",ylab="kappa")
abline(0,1,lty=2,col=2)
plot(NA,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",axes=FALSE)
}
