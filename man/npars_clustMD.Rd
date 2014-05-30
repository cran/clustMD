\name{npars_clustMD}
\alias{npars_clustMD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Calculates the number of free parameters for the clustMD model.
}
\description{
	Internal function.
}
\usage{
npars_clustMD(model, D, G, J, OrdIndx)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
The clustMD model fitted.
}
  \item{D}{
The dimension of the latent data.
}
  \item{G}{
The number of mixture components.
}
  \item{J}{
The number of variables.
}
  \item{OrdIndx}{
The number of continuous and ordinal (including binary) variables.
}
}
\details{
npars_clustMD: an internal function.
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