% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustMD_UtilityFunctions.R
\name{stable.probs}
\alias{stable.probs}
\title{Stable computation of the log of a sum}
\usage{
stable.probs(s)
}
\arguments{
\item{s}{a numeric vector.}
}
\value{
The log of the sum of the elements of \code{s}
}
\description{
Function takes a numeric vector and returns the log of the sum of the
elements of that vector. Calculations are done on the log scale for 
stability.
}
\keyword{internal}
