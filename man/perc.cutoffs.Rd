\name{perc.cutoffs}
\alias{perc.cutoffs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Calculates the threshold parameters for ordinal variables.
}
\description{
Internal function.
}
\usage{
perc.cutoffs(CnsIndx, OrdIndx, Y, N)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CnsIndx}{
The number of continuous variables.
}
  \item{OrdIndx}{
The number of continuous and ordinal (including binary) variables.
}
  \item{Y}{
An N x J data matrix.
}
  \item{N}{
Number of observations.
}
}
\details{
perc.cutoffs: an internal function.
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