% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustMD_InternalFunctions.R
\name{ObsLogLikelihood}
\alias{ObsLogLikelihood}
\title{Approximates the observed log likelihood.}
\usage{
ObsLogLikelihood(N, CnsIndx, G, Y, mu, Sigma, pi.vec, patt.indx, zlimits, J,
  OrdIndx, probs.nom, model, perc.cut, K)
}
\arguments{
\item{N}{the number of observations.}

\item{CnsIndx}{the number of continuous variables.}

\item{G}{the number of mixture components.}

\item{Y}{an \code{N x J} data matrix.}

\item{mu}{a \code{D x G} matrix of means.}

\item{Sigma}{a \code{D x D x G} array of covariance parameters.}

\item{pi.vec}{the mixing weights.}

\item{patt.indx}{a list of length equal to the number of observed response 
patterns. Each entry of the list details the observations for which that
response pattern was observed.}

\item{zlimits}{the truncation points for the latent data.}

\item{J}{the number of variables.}

\item{OrdIndx}{the sum of the number of continuous and ordinal (including
binary) variables.}

\item{probs.nom}{an array containing the response probabilities for each
nominal variable for each cluster}

\item{model}{the covariance model fitted to the data.}

\item{perc.cut}{threshold parameters.}

\item{K}{the number of levels for each variable.}
}
\value{
Output required for \code{clustMD} function.
}
\description{
Approximates the observed log likelihood.
}
\seealso{
\code{\link{clustMD}}
}
\keyword{internal}
