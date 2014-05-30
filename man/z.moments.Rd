\name{z.moments}
\alias{z.moments}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Calculates the first and second moments of the latent data.
}
\description{
	Internal function.
}
\usage{
z.moments(D, G, N, CnsIndx, OrdIndx, zlimits, mu, Sigma, Y, J, K, norms, nom.ind.Z)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{D}{
Dimension of the latent data.
}
  \item{G}{
	Number of mixture components.
}
  \item{N}{
Number of observations.
}
  \item{CnsIndx}{
The number of continuous variables.
}
  \item{OrdIndx}{
The number of continuous and ordinal (including binary) variables.
}
  \item{zlimits}{
	The truncation points for the latent data.
}
  \item{mu}{
A D x G matrix of means.
}
  \item{Sigma}{
A D x D x G array of covariance parameters.
}
  \item{Y}{
An N x J data matrix.
}
  \item{J}{
The number of variables.
}
  \item{K}{
	The number of levels for each variable.
}
  \item{norms}{
	A matrix of standard normal deviates.
}
  \item{nom.ind.Z}{
	The latent dimensions corresponding to each nominal variable.
}
}
\details{
z.moments: an internal function.
}
\value{
	Output required for clustMD function.
}
\references{
	McParland, D. and Gormley, I.C. (2014). Model based clustering for mixed data: clustMD. Technical report, University College Dublin.
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

