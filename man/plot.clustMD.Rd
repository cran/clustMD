% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustMD_Methods.R
\name{plot.clustMD}
\alias{plot.clustMD}
\title{Plotting method for objects of class \code{clustMD}}
\usage{
\method{plot}{clustMD}(x, ...)
}
\arguments{
\item{x}{a \code{clustMD} object.}

\item{...}{further arguments passed to or from other methods.}
}
\value{
Prints graphical summaries of the fitted model as detailed above.
}
\description{
Plots a parallel coordinates plot and dot plot of the estimated cluster 
means, a barplot of the variances by cluster for diagonal covariance models
or a heatmap of the covariance matrix for non-diagonal covariance 
structures, and a histogram of the clustering uncertainties for each 
observation.
}
\references{
McParland, D. and Gormley, I.C. (2016). Model based clustering 
    for mixed data: clustMD. Advances in Data Analysis and Classification, 
    10 (2):155-169.
}
\seealso{
\code{\link{clustMD}}
}
\keyword{device}
