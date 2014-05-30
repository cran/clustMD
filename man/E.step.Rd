\name{E.step}
\alias{E.step}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	E-step of the (MC)EM algorithm.
}
\description{
	Internal function.
}
\usage{
E.step(N, G, pi.vec, Y, OrdIndx, CnsIndx, D, perc.cut, mu, Sigma, Ez, J, probs.nom, K)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
	Number of observations.
}
  \item{G}{
	Number of mixture components.
}
  \item{pi.vec}{
	Mixing weights.
}
  \item{Y}{
	An N x J data matrix.
}
  \item{OrdIndx}{
	The number of continuous and ordinal (including binary) variables.
}
  \item{CnsIndx}{
	The number of continuous variables.
}
  \item{D}{
	Dimension of the latent data.
}
  \item{perc.cut}{
	Threshold parameters.
}
  \item{mu}{
	A D x G matrix of means.
}
  \item{Sigma}{
	A D x D x G array of covariance parameters.
}
  \item{Ez}{
	A N x D x G array of latent expected values.
}
  \item{J}{
	The number of observed variables.
}
  \item{probs.nom}{
	Probabilities of nominal responses by group.
}
  \item{K}{
	The number of levels for each variable.
}
}
\details{
	E-step: an internal function.
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

