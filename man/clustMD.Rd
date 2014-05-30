\name{clustMD}
\alias{clustMD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Model-Based Clustering for Mixed Data
}
\description{
	A function which fits the clustMD model to a data set consisting of any combination of continuous, binary, ordinal and nominal variables.
}
\usage{
clustMD(X, G, CnsIndx, OrdIndx, Nnorms, MaxIter, model, store.params = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
	A data matrix where the variables are ordered so that the continuous variables come first, the binary (coded 1 and 2) and ordinal variables (coded 1, 2,...) come second and the nominal variables (coded 1, 2,...) are in last position. 
}
  \item{G}{
	The number of mixture components to be fitted.
}
  \item{CnsIndx}{
	The number of continuous variables in the data set.
}
  \item{OrdIndx}{
	The sum of the number of continuous, binary and ordinal variables in the data set.
}
  \item{Nnorms}{
	The number of Monte Carlo samples to be used for the intractable E-step in the presence of nominal data.
}
  \item{MaxIter}{
	The number of iterations for which the (MC)EM algorithm should run.
}
  \item{model}{
	A string indicating which clustMD model is to be fitted. This may be one of: EII, VII, EEI, VEI, EVI or VVI.
}
  \item{store.params}{
	A logical variable indicating if the parameter estimates at each iteration should be saved and returned by the clustMD function.
}
}
\details{
Model-based clustering of mixed data using a parsimonious mixture of latent Gaussian variables.
}
\value{
A list is returned:
  \item{cl }{The cluster to which each observation belongs.}
  \item{tau }{A N x G matrix of the conditional probabilities of each observation blonging to each cluster.}
  \item{means }{A D x G matrix of the cluster means.}
  \item{A }{A G x D matrix containing the diagonal entries of the A matrix corresponding to each cluster.}
  \item{Lambda }{A G x D matrix of volume parameters corresponding to each observed or latent dimension for each cluster.}
  \item{Sigma}{A D x D x G array of the covariance matrices for each cluster.}
  \item{BIChat}{The estimated Bayesian information criterion for the model fitted.}
  \item{paramlist}{If store.params is true then paramlist is a list of the stored parameter values in the order given above with the saved estimated likelihood values in last position.}
}
\references{
	McParland, D. and Gormley, I.C. (2014). Model based clustering for mixed data: clustMD. Technical report, University College Dublin.
}
\author{
	Damien McParland
}

\examples{
	data(Byar)
	
	# Transformation skewed variables
Byar$Size.of.primary.tumour <- sqrt(Byar$Size.of.primary.tumour)
Byar$Serum.prostatic.acid.phosphatase <- log(Byar$Serum.prostatic.acid.phosphatase)

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

\dontrun{
	res <- clustMD(X=Y, G=3, CnsIndx=8, OrdIndx=11, Nnorms=20000, 
	MaxIter=100, model="EVI", store.params=FALSE)
}

}