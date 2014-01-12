\name{multiscale.rips}
\alias{multiscale.rips.ipca}
\title{Multiscale Rips Persistent Homology}
\description{
A multiscale approach to compute Rips persistent homology.
}


\usage{
multiscale.rips.ipca( X, epsGMRA,  maxD, single = FALSE )
}

\arguments{
  \item{X}{Data set. Each row is a data point}
  \item{epsGMRA}{Accuracy of GMRA approximation}
  \item{maxD}{Highest dimensionality of simplicies to insert into the filtration}
  \item{single}{Compute single scale TRUE or FALSE} 

}

\value{
  Persistent diagram, matrix with first column birth, second column death times
  and a third colum at which scale the death happened.
}

\author{
  Samuel Gerber
}

\examples{

##Noisy circle example
phi <- runif(100)*2*pi
arc <- cbind(cos(phi), sin(phi)) * (1+rnorm(length(phi)) * 0.1)

dgm <- multiscale.rips.ipca(X = arc, epsGMRA = 0, maxD = 2);

plot(x = dgm[,1], y = dgm[,2], xlim = c(0, max(dgm[,1:2]) ), ylim = c(0, max(dgm[,1:2])) ,
    asp=1, pch=19)
abline(0, 1)

}

\keyword{multiscale,persistence,homology,topology}
