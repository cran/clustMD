\name{M.step}
\alias{M.step}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	M-step of the (MC)EM algorithm.}
\description{
	Internal function.
}
\usage{
M.step(tau, N, Elz, J, OrdIndx, D, G, Y, CnsIndx, S2, model, a)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tau}{
	An N x G matrix of conditional cluster membership probabilities.
}
  \item{N}{
	Number of observations.	
}
  \item{Elz}{
	A N x D x G array of expected values of cluster labels multiplied by the latent data.
}
  \item{J}{
	The number of variables.
}
  \item{OrdIndx}{
	The number of continuous and ordinal (including binary) variables.}
  \item{D}{
	Dimension of the latent data.
}
  \item{G}{
	The number of mixture components.
}
  \item{Y}{
	A N x J data matrix.
}
  \item{CnsIndx}{
	The number of continuous variables.
}
  \item{S2}{
	A N x D x G array of expected values of the latent data squared.
}
  \item{model}{
	Which clustMD model is fitted.
}
  \item{a}{
	A G x D matrix of the entries of A.
}
}
\details{
	M-step: an internal function.}
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