% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustMD_InternalFunctions.R
\name{E.step}
\alias{E.step}
\title{E-step of the (MC)EM algorithm}
\usage{
E.step(N, G, D, CnsIndx, OrdIndx, zlimits, mu, Sigma, Y, J, K, norms, nom.ind.Z,
  patt.indx, pi.vec, model, perc.cut)
}
\arguments{
\item{N}{number of observations.}

\item{G}{number of mixture components.}

\item{D}{dimension of the latent data.}

\item{CnsIndx}{the number of continuous variables.}

\item{OrdIndx}{the sum of the number of continuous and ordinal (including
binary) variables.}

\item{zlimits}{the truncation points for the latent data.}

\item{mu}{a D x G matrix of means.}

\item{Sigma}{a D x D x G array of covariance parameters.}

\item{Y}{an N x J data matrix.}

\item{J}{the number of observed variables.}

\item{K}{the number of levels for each variable.}

\item{norms}{a matrix of standard normal deviates.}

\item{nom.ind.Z}{the latent dimensions corresponding to each nominal 
variable.}

\item{patt.indx}{a list of length equal to the number of observed response 
patterns. Each entry of the list details the observations for which that
response pattern was observed.}

\item{pi.vec}{mixing weights.}

\item{model}{the covariance model fitted to the data.}

\item{perc.cut}{threshold parameters.}
}
\value{
Output required for clustMD function.
}
\description{
Internal function.
}
\seealso{
\code{\link{clustMD}}
}
\keyword{internal}
