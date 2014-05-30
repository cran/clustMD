\name{ObsLogLikelihood}
\alias{ObsLogLikelihood}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Approximates the observed log likelihood.
}
\description{
	Internal function.
}
\usage{
ObsLogLikelihood(N, CnsIndx, G, Y, mu, Sigma, pi.vec, J, 
OrdIndx, K, perc.cut, Nnorms, zlimits, nom.ind.Z)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
	The number of observations.
}
  \item{CnsIndx}{
	The number of continuous variables.
}
  \item{G}{
	The number of mixture components.
}
  \item{Y}{
	An N x J data matrix.
}
  \item{mu}{
	A D x G matrix of means.
}
  \item{Sigma}{
	A D x D x G array of covariance parameters.
}
  \item{pi.vec}{
	The mixing weights.
}
  \item{J}{
The number of variables.
}
  \item{OrdIndx}{
The number of continuous and ordinal (including binary) variables.
}
  \item{K}{
The number of levels for each variable.
}
  \item{perc.cut}{
Threshold parameters.
}
  \item{Nnorms}{
The number of Monte Carlo samples.
}
  \item{zlimits}{
The truncation points given each response.
}
  \item{nom.ind.Z}{
A list indicating the latent dimensions corresponding to each nominal variable.
}
}
\details{
	ObsLogLikelihood: an internal function.
}
\value{
Output required for clustMD function.
}
\references{
	McParland, D. and Gormley, I.C. (2014). Model based clustering for mixed data: 
	clustMD. Technical report, University College Dublin.
}
\author{
	Damien McParland
}
\note{
	An internal function.
}

\seealso{
	\code{\link{clustMD}}
}

\keyword{ internal }

