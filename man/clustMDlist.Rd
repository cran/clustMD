% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustMD_UserFunctions.R
\name{clustMDlist}
\alias{clustMDlist}
\title{Model Based Clustering for Mixed Data}
\usage{
clustMDlist(arglist)
}
\arguments{
\item{arglist}{a list of input arguments for \code{clustMD}. See 
\code{\link{clustMD}}.}
}
\value{
A \code{clustMD} object. See \code{\link{clustMD}}.
}
\description{
A function that fits the clustMD model to a data set consisting of any
combination of continuous, binary, ordinal and nominal variables. This
function is a wrapper for \code{\link{clustMD}} that takes arguments as a
list.
}
\examples{
    data(Byar)

    # Transformation skewed variables
    Byar$Size.of.primary.tumour <- sqrt(Byar$Size.of.primary.tumour)
    Byar$Serum.prostatic.acid.phosphatase <- 
    log(Byar$Serum.prostatic.acid.phosphatase)

    # Order variables (Continuous, ordinal, nominal)
    Y <- as.matrix(Byar[, c(1, 2, 5, 6, 8, 9, 10, 11, 3, 4, 12, 7)])

    # Start categorical variables at 1 rather than 0
    Y[, 9:12] <- Y[, 9:12] + 1

    # Standardise continuous variables
    Y[, 1:8] <- scale(Y[, 1:8])

    # Merge categories of EKG variable for efficiency
    Yekg <- rep(NA, nrow(Y))
    Yekg[Y[,12]==1] <- 1
    Yekg[(Y[,12]==2)|(Y[,12]==3)|(Y[,12]==4)] <- 2
    Yekg[(Y[,12]==5)|(Y[,12]==6)|(Y[,12]==7)] <- 3
    Y[, 12] <- Yekg

    argList <- list(X=Y, G=3, CnsIndx=8, OrdIndx=11, Nnorms=20000,
    MaxIter=500, model="EVI", store.params=FALSE, scale=TRUE, 
    startCL="kmeans", autoStop=FALSE, ma.band=50, stop.tol=NA)

    \dontrun{
    res <- clustMDlist(argList)
    }

}
\references{
McParland, D. and Gormley, I.C. (2016). Model based clustering
    for mixed data: clustMD. Advances in Data Analysis and Classification,
    10 (2):155-169.
}
\seealso{
\code{\link{clustMD}}
}
