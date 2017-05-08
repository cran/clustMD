# This file contains functions that will be called by the end user of clustMD.

#' Model Based Clustering for Mixed Data
#'
#' A function that fits the clustMD model to a data set consisting of any
#' combination of continuous, binary, ordinal and nominal variables.
#'
#' @param X a data matrix where the variables are ordered so that the 
#'     continuous variables come first, the binary (coded 1 and 2) and
#'     ordinal variables (coded 1, 2, ...) come second and the nominal
#'     variables (coded 1, 2, ...) are in last position.
#' @param G the number of mixture components to be fitted.
#' @param CnsIndx the number of continuous variables in the data set.
#' @param OrdIndx the sum of the number of continuous, binary and ordinal
#'     variables in the data set.
#' @param Nnorms the number of Monte Carlo samples to be used for the 
#'     intractable E-step in the presence of nominal data. Irrelevant if
#'     there are no nominal variables.
#' @param MaxIter the maximum number of iterations for which the (MC)EM 
#'     algorithm should run.
#' @param model a string indicating which clustMD model is to be fitted. This
#'     may be one of: \code{EII, VII, EEI, VEI, EVI, VVI} or \code{BD}.
#' @param store.params a logical argument indicating if the parameter estimates
#'     at each iteration should be saved and returned by the clustMD function.
#' @param scale a logical argument indicating if the continuous variables
#'     should be standardised.
#' @param startCL a string indicating which clustering method should be used to
#'     initialise the (MC)EM algorithm. This may be one of "kmeans" (K means
#'     clustering), "hclust" (hierarchical clustering), "mclust" (finite 
#'     mixture of Gaussian distributions), "hc_mclust" (model-based 
#'     hierarchical clustering) or "random" (random cluster allocation).
#' @param autoStop a logical argument indicating whether the (MC)EM algorithm
#'     should use a stopping criterion to decide if convergence has been 
#'     reached. Otherwise the algorithm will run for \code{MaxIter} iterations. 
#' 
#'     If only continuous variables are present the algorithm will use Aitken's
#'     acceleration criterion with tolerance \code{stop.tol}.
#' 
#'     If categorical variables are present, the stopping criterion is based
#'     on a moving average of the approximated log likelihood values. Let 
#'     \code{t} denote the current interation. The average of the 
#'     \code{ma.band} most recent approximated log likelihood values is 
#'     compared to the average of another \code{ma.band} iterations with a lag
#'     of 10 iterations. If this difference is less than the tolerance the
#'     algorithm will be said to have converged.
#' @param ma.band the number of iterations to be included in the moving average
#'     calculation for the stopping criterion. 
#' @param stop.tol the tolerance of the (MC)EM stopping criterion.
#'
#' @return An object of class clustMD is returned. The output components are as
#'     follows:
#'     \item{model }{The covariance model fitted to the data.}
#'     \item{G }{The number of clusters fitted to the data.}
#'     \item{Y }{The observed data matrix.}
#'     \item{cl }{The cluster to which each observation belongs.}
#'     \item{tau }{A \code{N x G} matrix of the probabilities of 
#'         each observation blonging to each cluster.}
#'     \item{means }{A \code{D x G} matrix of the cluster means. Where D is the 
#'         dimension of the combined observed and latent continuous space.}
#'     \item{A }{A \code{G x D} matrix containing the diagonal entries of the 
#'         \code{A} matrix corresponding to each cluster.}
#'     \item{Lambda }{A \code{G x D} matrix of volume parameters corresponding
#'         to each observed or latent dimension for each cluster.}
#'     \item{Sigma }{A \code{D x D x G} array of the covariance matrices for 
#'         each cluster.}
#'     \item{BIChat }{The estimated Bayesian information criterion for the 
#'         model fitted.}
#'     \item{ICLhat }{The estimated integrated classification likelihood criterion
#'         for the model fitted.}
#'     \item{paramlist }{If store.params is \code{TRUE} then paramlist is a 
#'         list of the stored parameter values in the order given above with 
#'         the saved estimated likelihood values in last position.}
#'     \item{Varnames }{A character vector of names corresponding to the 
#'         columns of \code{Y}}
#'     \item{Varnames_sht }{A truncated version of \code{Varnames}. Used for
#'         plotting.}
#'     \item{likelihood.store }{A vector containing the estimated log 
#'         likelihood at each iteration.}
#' 
#' @references McParland, D. and Gormley, I.C. (2016). Model based clustering
#'     for mixed data: clustMD. Advances in Data Analysis and Classification,
#'     10 (2):155-169.
#' @export
#' 
#' @importFrom mclust hcVVV
#' @importFrom mclust mclustBIC
#' @import utils
#' @import stats
#' 
#' @examples 	data(Byar)
#'     # Transformation skewed variables
#'     Byar$Size.of.primary.tumour <- sqrt(Byar$Size.of.primary.tumour)
#'     Byar$Serum.prostatic.acid.phosphatase <- log(Byar$Serum.prostatic.acid.phosphatase)
#' 
#'     # Order variables (Continuous, ordinal, nominal)
#'     Y <- as.matrix(Byar[, c(1, 2, 5, 6, 8, 9, 10, 11, 3, 4, 12, 7)])
#' 
#'     # Start categorical variables at 1 rather than 0
#'     Y[, 9:12] <- Y[, 9:12] + 1
#' 
#'     # Standardise continuous variables
#'     Y[, 1:8] <- scale(Y[, 1:8])
#' 
#'     # Merge categories of EKG variable for efficiency
#'     Yekg <- rep(NA, nrow(Y))
#'     Yekg[Y[,12]==1] <- 1
#'     Yekg[(Y[,12]==2)|(Y[,12]==3)|(Y[,12]==4)] <- 2
#'     Yekg[(Y[,12]==5)|(Y[,12]==6)|(Y[,12]==7)] <- 3
#'     Y[, 12] <- Yekg
#' 
#'     \dontrun{
#'     res <- clustMD(X = Y, G = 3, CnsIndx = 8, OrdIndx = 11, Nnorms = 20000,
#'     MaxIter = 500, model = "EVI", store.params = FALSE, scale = TRUE, 
#'     startCL = "kmeans", autoStop= TRUE, ma.band=30, stop.tol=0.0001)
#'     }

clustMD <- function(X, G, CnsIndx, OrdIndx, Nnorms, MaxIter, model, 
                    store.params = FALSE, scale = FALSE, startCL = "hc_mclust",
                    autoStop=FALSE, ma.band=50, stop.tol=NA) {
  # Controls
  Y <- as.matrix(X)
  N <- nrow(Y)
  J <- ncol(Y)
  
  # Scale continuous data if required
  if (scale) {
    if (CnsIndx > 0) 
      Y[, 1:CnsIndx] <- scale(Y[, 1:CnsIndx])
  }  # if
  
  # Number of levels for each categorical item
  K <- apply(Y, 2, max)
  if (CnsIndx > 0) 
    K[1:CnsIndx] <- NA
  
  # Dimension of latent space
  D <- J
  if (J > OrdIndx) 
    D <- OrdIndx + sum(K[(OrdIndx + 1):J] - 1)
  
  # Which dimensions correspond to each item
  if (J > OrdIndx) {
    nom.ind.Z <- vector("list", J - OrdIndx)
    for (j in 1:(J - OrdIndx)) {
      if (j == 1) {
        start <- OrdIndx + 1
      } else {
        start <- OrdIndx + sum(K[(OrdIndx + 1):(OrdIndx + j - 1)] - 1) + 1
      }
      finish <- start + K[OrdIndx + j] - 2
      nom.ind.Z[[j]] <- c(start:finish)
    }  # j
  }  # if
  
  if ( (model == "BD")&(OrdIndx > CnsIndx) ){
    if((OrdIndx-CnsIndx)==1){
      patt.indx <- list()
      for(p in 1:max(Y[, OrdIndx])) {patt.indx[[p]] <- which(Y[, OrdIndx]==p)}
    }else{
      patt.tab <- data.frame(table(data.frame((Y[, (CnsIndx + 1):OrdIndx]))))
      patt.tab <- patt.tab[patt.tab$Freq != 0, 1:(OrdIndx - CnsIndx)]
      patt.indx <- list()
      for (p in 1:nrow(patt.tab)) {
        patt.indx[[p]] <- which(apply(Y[, (CnsIndx + 1):OrdIndx], 1, patt.equal, patt.tab[p, ]))
      } # p
    }
  } # if
  
  
  # Names for the observed/latent continuous dimensions
  # Long version
  Ynames <- colnames(X)
  VarNames <- as.character(1:D)
  if (!is.null(Ynames)) {
    VarNames[1:OrdIndx] <- Ynames[1:OrdIndx]
    
    if (J > OrdIndx) {
      NomNames <- list()
      for (j in (OrdIndx + 1):J) {
        NomNames[[j - OrdIndx]] <- rep(NA, (K[j] - 1))
        for(k in 1:(K[j] - 1)) {
          NomNames[[j - OrdIndx]][k] <- paste(Ynames[j],"_", k, sep = "")
        }
      }
      VarNames[(OrdIndx + 1):D] <- unlist(NomNames)
    }
  } else {
    for(j in 1:OrdIndx) {
      VarNames[j] <- paste("V", j, sep = "")
    }
    
    if (J > OrdIndx) {
      NomNames <- list()
      for (j in (OrdIndx + 1):J) {
        NomNames[[j - OrdIndx]] <- rep(paste("V", j, sep=""), (K[j] - 1))
        for(k in 1:(K[j] - 1)) {
          NomNames[[j - OrdIndx]][k] <- paste("V", j, "_", k, sep = "")
        }
      }
      VarNames[(OrdIndx + 1):D] <- unlist(NomNames)
    }
  }
  
  # Shortened version
  VarNames_sht <- as.character(1:D)
  if (!is.null(Ynames)) {
    VarNames_sht[1:OrdIndx] <- substr(Ynames[1:OrdIndx], 1, 7)
    
    if (J > OrdIndx) {
      NomNames_sht <- list()
      for (j in (OrdIndx + 1):J) {
        NomNames_sht[[j - OrdIndx]] <- rep(NA, (K[j] - 1))
        for(k in 1:(K[j] - 1)) {
          NomNames_sht[[j - OrdIndx]][k] <- paste(substr(Ynames[j], 1, 7),"_", k, sep = "")
        }
      }
      VarNames_sht[(OrdIndx + 1):D] <- unlist(NomNames_sht)
    }
  } else {
    VarNames_sht <- VarNames
  }
  
  
  ### Initial Values
  
  ### Estimated starting values Expected value of latent data
  Ez <- array(NA, c(N, D, G))
  for (g in 1:G) Ez[, 1:J, g] <- Y
  
  # Cutoffs for ordinal items
  if (OrdIndx > CnsIndx) {
    perc.cut <- perc.cutoffs(CnsIndx, OrdIndx, Y, N)
    zlimits <- array(NA, c(N, J, 2))
    zlimits[, 1:CnsIndx, 1] <- -Inf
    zlimits[, 1:CnsIndx, 2] <- Inf
    for (j in (CnsIndx + 1):OrdIndx) {
      for (k in 1:K[j]) {
        zlimits[Y[, j] == k, j, 1] <- perc.cut[[j]][k]
        zlimits[Y[, j] == k, j, 2] <- perc.cut[[j]][k + 1]
      }
    }
  } else {
    perc.cut <- list()
    zlimits <- array(NA, c(N, J, 2))
  }
  
  # Initial latent data matrix
  Zstart <- function(Kj, y) {
    new.z <- rep(0, (Kj - 1))
    if (y == 1) {
      new.z <- msm::rtnorm((Kj - 1), mean = 0, sd = 1, upper = 0)
    } else {
      new.z[-(y - 1)] <- stats::rnorm((Kj - 2), mean = 0, sd = 1)
      new.z[(y - 1)] <- msm::rtnorm(1, mean = 0, sd = 1, lower = max(new.z))
    }
    new.z
  }
  
  Zinit <- matrix(NA, N, D)
  Zinit[, 1:OrdIndx] <- Y[, 1:OrdIndx]
  if (J > OrdIndx) {
    for (j in (OrdIndx + 1):J) {
      for (i in 1:N) {
        Zinit[i, nom.ind.Z[[j - OrdIndx]]] <- Zstart(K[j], Y[i, j])
      }  # i
    }  # j
  }
  
  # Initial clustering
  if (startCL == "kmeans") {
    if (CnsIndx > 0) {
      ind <- kmeans(Y[, 1:CnsIndx], G)$cl
    } else {
      ind <- kmeans(Y, G)$cl
    }  # if
    
  } else if (startCL == "hclust") {
    temp <- hclust(dist(Y))
    ind <- cutree(temp, G)
    
  } else if (startCL == "mclust") {
    if (CnsIndx > 0) {
      if(CnsIndx==1){
        ind <- mclust::Mclust(Y[, 1:CnsIndx], G, "V")$cl
      } else {
        ind <- mclust::Mclust(Y[, 1:CnsIndx], G, "VVV")$cl
      } # ifelse
    } else {
      if(J==1){
        ind <- mclust::Mclust(Y, G, "V")$cl
      } else {
        ind <- mclust::Mclust(Y, G, "VVV")$cl
      } # ifelse
    }
    
  } else if (startCL == "random") {
    ind <- sample(1:G, N, replace = TRUE)
    
  } else if (startCL == "hc_mclust") {
    if (CnsIndx > 0) {
      ind <- mclust::hclass(mclust::hc(modelName = "VVV", data = Y[, 1:CnsIndx]), G)
    } else {
      stop("Cannot use hc_mclust since there are no continuous variables")
    }
  } else {
    stop("Unknown starting value algorithm chosen! \n 
        Choose from: kmeans, hclust, hc_mclust, mclust or random")
  }
  
  # Mixing weights
  pi.vec <- table(ind)/N
  
  # Mean
  mu <- matrix(NA, D, G)
  for (g in 1:G) 
    mu[, g] <- colMeans(matrix(Zinit[ind == g, ], sum(ind == g), D))
  
  # Covaraince
  Sigma <- array(NA, c(D, D, G))
  for (g in 1:G) Sigma[, , g] <- diag(D)
  
  # A matrix
  a <- matrix(1, G, D)
  
  likeStore <- rep(NA, MaxIter)
  
  ## Storage
  if (store.params == TRUE) {
    ind.store <- matrix(NA, N, MaxIter)
    Ez.store <- array(NA, c(N, D, G, MaxIter))
    tau.store <- array(NA, c(N, G, MaxIter))
    mu.store <- array(NA, c(D, G, MaxIter))
    lambda.store <- array(NA, c(G, D, MaxIter))
    a.store <- array(NA, c(G, D, MaxIter))
    Sigma.store <- array(NA, c(D, D, G, MaxIter))
    # likeStore <- rep(NA, MaxIter)
    if (J > OrdIndx) 
      probs.nom.store <- array(NA, c(J - OrdIndx, max(K[(OrdIndx + 1):J]), 
                                     G, MaxIter))
  }
  
  pb <- txtProgressBar(style = 3)
  prog <- 1
  ### EM Loop
  for (iter in 1:MaxIter) {
    # if(iter%%10==0) print(iter)
    if (iter > prog * MaxIter/10) 
      prog <- prog + 1
    setTxtProgressBar(pb, prog * 0.1)
    
    # Standard normal deviates for MC approximation
    if (J > OrdIndx) 
      norms <- MASS::mvrnorm(Nnorms, mu = rep(0, max(K[(OrdIndx + 1):J]) - 1),
                             Sigma = diag(max(K[(OrdIndx + 1):J]) - 1))
    
    # E-step
    temp.E <- E.step(N, G, D, CnsIndx, OrdIndx, zlimits, mu, Sigma, Y, 
                     J, K, norms, nom.ind.Z, patt.indx, pi.vec, model, perc.cut)
    
    tau <- temp.E[[1]]
    sumTauEz <- temp.E[[2]]
    sumTauS <- temp.E[[3]]
    probs.nom <- temp.E[[4]]
    Ez <- temp.E[[5]]   # ICL test
    ind <- mclust::map(tau)
    
    # M-step
    temp.M <- M.step(tau, N, sumTauEz, J, OrdIndx, D, G, Y, CnsIndx, sumTauS, 
                     model, a, nom.ind.Z)
    pi.vec <- temp.M[[1]]
    mu <- temp.M[[2]]
    lambda <- temp.M[[3]]
    a <- temp.M[[4]]
    Sigma <- temp.M[[5]]
    
    likeStore[iter] <- ObsLogLikelihood(N, CnsIndx, G, Y, mu, Sigma, pi.vec,
                                        patt.indx, zlimits, J, OrdIndx, probs.nom,
                                        model, perc.cut, K)
    
    # Store
    if (store.params == TRUE) {
      ind.store[, iter] <- ind
      # Ez.store[, , ,iter] <- Ez
      tau.store[, , iter] <- tau
      mu.store[, , iter] <- mu
      lambda.store[, , iter] <- lambda
      a.store[, , iter] <- a
      Sigma.store[, , , iter] <- Sigma
    }
    
    ### Check stopping criterion
    if(autoStop & (iter > 5)) { 
      autoStop.aitken <- (OrdIndx==J)
      autoStop.ma <- (OrdIndx < J) & (iter > (ma.band+10)) & (iter%%10==0)
      
      # Only continuous and ordinal present
      if (autoStop.aitken) {
        # Check input
        if(is.na(stop.tol)) { stop("stop.tol not specified.") }
        
        check.diff <- (likeStore[iter] - likeStore[iter-1])/abs(1 + likeStore[iter])
        if( check.diff < stop.tol) {
          setTxtProgressBar(pb, 1)
          break
        }
      } # if
    
      # Nominal variable present
      if (autoStop.ma) {
        # Check inputs
        if(is.na(stop.tol)) { stop("stop.tol or ma.band not specified.") }
        
        MA.lag <- 10
        likeOld <- mean(likeStore[(iter - ma.band - MA.lag):(iter - MA.lag)])
        likeNew <- mean(likeStore[(iter - ma.band):iter])
        
        if(abs((likeNew - likeOld)/likeNew) < stop.tol) { 
          setTxtProgressBar(pb, 1)
          break 
        }  # if
      } # if
    } # if
  }  # iter
  
  # approximated BIC
  if (model == "BD"){
    probs.nom <- z.moments(D, G, N, CnsIndx, OrdIndx, zlimits, mu, Sigma, 
                           Y, J, K, norms, nom.ind.Z, patt.indx)[[2]]
  } else {
    probs.nom <- z.moments_diag(D, G, N, CnsIndx, OrdIndx, zlimits, mu, Sigma, 
                                Y, J, K, norms, nom.ind.Z)[[2]]
  }
  
  obslike <- ObsLogLikelihood(N, CnsIndx, G, Y, mu, Sigma, pi.vec, patt.indx, 
                              zlimits, J, OrdIndx, probs.nom, model, perc.cut, K)
  
  BIChat <- 2 * obslike - npars_clustMD(model, D, G, J, CnsIndx, OrdIndx, K) * log(N)
  close(pb)
  
  # approximated ICL
  CompleteLike_i <- rep(NA, N)
  for(i in 1:N){
    CompleteLike_i[i] <- log(pi.vec[ind[i]]) +
      mvtnorm::dmvnorm(Ez[i, , ind[i]], mean=mu[, ind[i]], 
                       sigma=matrix(Sigma[, , ind[i]], nrow = D, ncol = D), log = TRUE)
  }
  
  ICLhat <- 2*sum(CompleteLike_i) - npars_clustMD(model, D, G, J, CnsIndx, OrdIndx, K) * log(N)
  
  
  ## Give variable names to output matrices mu
  rownames(mu) <- VarNames
  # Sigma
  rownames(Sigma) <- VarNames
  colnames(Sigma) <- VarNames
  
  if (store.params == TRUE) {
    params.store.list <- list(cl.store = ind.store, tau.store = tau.store, 
                              means.store = mu.store, A.store = a.store, 
                              lambda.store = lambda.store, Sigma.store = Sigma.store)
  } else {
    params.store.list <- NULL
  }
  out.clustMD <- list(model = model, G = G, Y=Y, cl = ind, tau = tau, means = mu, 
                      A = a, Lambda = lambda, Sigma = Sigma, BIChat = BIChat,
                      ICLhat=ICLhat, paramlist = params.store.list, VarNames=VarNames,
                      VarNames_sht=VarNames_sht, likelihood.store = likeStore)
  class(out.clustMD) <- "clustMD"
  out.clustMD
}

#-------------------------------------------------------------------------#

### Wrapper function for clustMD that takes inputs for several models as a
### list Will use in clustMDparallel

#' Model Based Clustering for Mixed Data
#' 
#' A function that fits the clustMD model to a data set consisting of any
#' combination of continuous, binary, ordinal and nominal variables. This
#' function is a wrapper for \code{\link{clustMD}} that takes arguments as a
#' list.
#'
#' @param arglist a list of input arguments for \code{clustMD}. See 
#'     \code{\link{clustMD}}.
#'
#' @return A \code{clustMD} object. See \code{\link{clustMD}}.
#'
#' @references McParland, D. and Gormley, I.C. (2016). Model based clustering
#'     for mixed data: clustMD. Advances in Data Analysis and Classification,
#'     10 (2):155-169.
#'
#' @examples 
#'     data(Byar)
#' 
#'     # Transformation skewed variables
#'     Byar$Size.of.primary.tumour <- sqrt(Byar$Size.of.primary.tumour)
#'     Byar$Serum.prostatic.acid.phosphatase <- 
#'     log(Byar$Serum.prostatic.acid.phosphatase)
#' 
#'     # Order variables (Continuous, ordinal, nominal)
#'     Y <- as.matrix(Byar[, c(1, 2, 5, 6, 8, 9, 10, 11, 3, 4, 12, 7)])
#' 
#'     # Start categorical variables at 1 rather than 0
#'     Y[, 9:12] <- Y[, 9:12] + 1
#' 
#'     # Standardise continuous variables
#'     Y[, 1:8] <- scale(Y[, 1:8])
#' 
#'     # Merge categories of EKG variable for efficiency
#'     Yekg <- rep(NA, nrow(Y))
#'     Yekg[Y[,12]==1] <- 1
#'     Yekg[(Y[,12]==2)|(Y[,12]==3)|(Y[,12]==4)] <- 2
#'     Yekg[(Y[,12]==5)|(Y[,12]==6)|(Y[,12]==7)] <- 3
#'     Y[, 12] <- Yekg
#' 
#'     argList <- list(X=Y, G=3, CnsIndx=8, OrdIndx=11, Nnorms=20000,
#'     MaxIter=500, model="EVI", store.params=FALSE, scale=TRUE, 
#'     startCL="kmeans", autoStop=FALSE, ma.band=50, stop.tol=NA)
#' 
#'     \dontrun{
#'     res <- clustMDlist(argList)
#'     }
#' 
#' @seealso \code{\link{clustMD}}
#' 
clustMDlist <- function(arglist) {
  ## Define arguments
  Y <- arglist[[1]]
  G <- arglist[[2]]
  CnsIndx <- arglist[[3]]
  OrdIndx <- arglist[[4]]
  Nnorms <- arglist[[5]]
  MaxIter <- arglist[[6]]
  model <- arglist[[7]]
  store.params <- arglist[[8]]
  scale <- arglist[[9]]
  startCL <- arglist[[10]]
  
  autoStop <- arglist[[11]]
  ma.band <- arglist[[12]]
  stop.tol <- arglist[[13]]
  
  # run clustMD
  # res <- tryCatch(clustMD(Y, G, CnsIndx, OrdIndx, Nnorms, MaxIter, model,
  #   store.params, scale, startCL), error = function(e) NULL)
  res <- tryCatch(clustMD(Y, G, CnsIndx, OrdIndx, Nnorms, MaxIter, model,
                          store.params, scale, startCL, autoStop, ma.band,
                          stop.tol), error = function(e) NULL)
  
  if (is.null(res)) {
    return(NULL)
  } else {
    return(res)
  }  # ifelse
}

#-------------------------------------------------------------------------#

### Function to fit multiple clustMD models in parallel Each model is fitted
### in parallel using teh parallel package. No parallelisation at lower levels than that.
### models is a character vector indicating the covariance models to be
### fitted G is a vector indicating the numbers of clusters to be fitted for
### each element of models

#' Run multiple clustMD models in parallel
#' 
#' This function allows the user to run multiple clustMD models in parallel.
#' The inputs are similar to \code{clustMD()} except \code{G} is now a vector
#' containing the the numbers of components the user would like to fit and 
#' \code{models} is a vector of strings indicating the covariance models the 
#' user would like to fit for each element of G. The user can specify the 
#' number of cores to be used or let the function detect the number available.
#'
#' @param X a data matrix where the variables are ordered so that the 
#'     continuous variables come first, the binary (coded 1 and 2) and ordinal
#'     variables (coded 1, 2,...) come second and the nominal variables 
#'     (coded 1, 2,...) are in last position.
#' @param CnsIndx the number of continuous variables in the data set.
#' @param OrdIndx the sum of the number of continuous, binary and ordinal 
#'     variables in the data set.
#' @param G a vector containing the numbers of mixture components to be fitted.
#' @param models a vector of strings indicating which clustMD models are to be
#'     fitted. This may be one of: \code{EII, VII, EEI, VEI, EVI, VVI} or 
#'     \code{BD}.
#' @param Nnorms the number of Monte Carlo samples to be used for the 
#'     intractable E-step in the presence of nominal data.
#' @param MaxIter the maximum number of iterations for which the (MC)EM 
#'     algorithm should run.
#' @param store.params a logical variable indicating if the parameter estimates
#'     at each iteration should be saved and returned by the \code{clustMD}
#'     function.
#' @param scale a logical variable indicating if the continuous variables 
#'     should be standardised.
#' @param startCL a string indicating which clustering method should be used to
#'     initialise the (MC)EM algorithm. This may be one of "kmeans" (K means 
#'     clustering), "hclust" (hierarchical clustering), "mclust" (finite 
#'     mixture of Gaussian distributions), "hc_mclust" (model-based 
#'     hierarchical clustering) or "random" (random cluster allocation).
#' @param Ncores the number of cores the user would like to use. Must be less
#'     than or equal to the number of cores available.
#' @param autoStop a logical argument indicating whether the (MC)EM algorithm
#'     should use a stopping criterion to decide if convergence has been 
#'     reached. Otherwise the algorithm will run for \code{MaxIter} iterations. 
#' 
#'     If only continuous variables are present the algorithm will use Aitken's
#'     acceleration criterion with tolerance \code{stop.tol}. 
#' 
#'     If categorical variables are present, the stopping criterion is based
#'     on a moving average of the approximated log likelihood values. let $t$
#'     denote the current interation. The average of the \code{ma.band} most
#'     recent approximated log likelihood values is compared to the average of
#'     another \code{ma.band} iterations with a lag of 10 iterations.
#'     If this difference is less than the tolerance the algorithm will be 
#'     said to have converged.
#' @param ma.band the number of iterations to be included in the moving average
#'     stopping criterion. 
#' @param stop.tol the tolerance of the (MC)EM stopping criterion.
#'
#' @return An object of class \code{clustMDparallel} is returned. The output 
#'     components are as follows:
#'     \item{BICarray }{A matrix indicating the estimated BIC values for each
#'         of the models fitted.}
#'     \item{results }{A list containing the output for each of the models 
#'         fitted. Each entry of this list is a \code{clustMD} object. If the 
#'         algorithm failed to fit a particular model, the corresponding entry 
#'         of \code{results} will be \code{NULL}.}
#' @export
#' 
#' @references McParland, D. and Gormley, I.C. (2016). Model based clustering 
#'     for mixed data: clustMD. Advances in Data Analysis and Classification, 
#'     10 (2):155-169.
#' 
#' @seealso \code{\link{clustMD}}
#'
#' @examples 
#'     data(Byar)
#' 
#'     # Transformation skewed variables
#'     Byar$Size.of.primary.tumour <- sqrt(Byar$Size.of.primary.tumour)
#'     Byar$Serum.prostatic.acid.phosphatase <- 
#'     log(Byar$Serum.prostatic.acid.phosphatase)
#' 
#'     # Order variables (Continuous, ordinal, nominal)
#'     Y <- as.matrix(Byar[, c(1, 2, 5, 6, 8, 9, 10, 11, 3, 4, 12, 7)])
#' 
#'     # Start categorical variables at 1 rather than 0
#'     Y[, 9:12] <- Y[, 9:12] + 1
#' 
#'     # Standardise continuous variables
#'     Y[, 1:8] <- scale(Y[, 1:8])
#' 
#'     # Merge categories of EKG variable for efficiency
#'     Yekg <- rep(NA, nrow(Y))
#'     Yekg[Y[,12]==1] <- 1
#'     Yekg[(Y[,12]==2)|(Y[,12]==3)|(Y[,12]==4)] <- 2
#'     Yekg[(Y[,12]==5)|(Y[,12]==6)|(Y[,12]==7)] <- 3
#'     Y[, 12] <- Yekg
#' 
#'     \dontrun{
#'     res <- clustMDparallel(X = Y, G = 1:3, CnsIndx = 8, OrdIndx = 11, Nnorms = 20000,
#'     MaxIter = 500, models = c("EVI", "EII", "VII"), store.params = FALSE, scale = TRUE, 
#'     startCL = "kmeans", autoStop= TRUE, ma.band=30, stop.tol=0.0001)
#'   
#'     res$BICarray
#' }
#' 
clustMDparallel <- function(X, CnsIndx, OrdIndx, G, models, Nnorms, MaxIter,
                            store.params, scale, startCL = "hc_mclust", Ncores = NULL,
                            autoStop = FALSE, ma.band = 50, stop.tol = NA) {
  
  ### Specify the number of cores to be in the 'cluster'
  Nmods <- length(G)*length(models)
  temp_cores <- parallel::detectCores()
  if (is.null(Ncores)) {
    cores <- min(temp_cores, Nmods)
  } else if (Ncores <= temp_cores) {
    cores <- min(Ncores, Nmods)
  } else {
    cat("Ncores is greater than the number available (", temp_cores, ") setting Ncores
        equal to ", min(temp_cores, Nmods))
    cores <- min(temp_cores, Nmods)
  }  # if
  
  ### Matrix to store estimated BIC of fitted models
  BICarray <- matrix(NA, length(G), length(models))
  rNames <- rep("", length(G))
  for (i in 1:length(G)) rNames[i] <- paste("G=", G[i], sep = "")
  rownames(BICarray) <- rNames
  colnames(BICarray) <- models
  
  ICLarray <- BICarray
  
  ### Create a list whose entries are themselves lists of inputs for clustMD
  ArgsList <- list()
  i <- 1
  for (m in 1:length(models)) {
    for (g in G) {
      # ArgsList[[i]] <- list(X, g, CnsIndx, OrdIndx, Nnorms, MaxIter, 
      #   models[m], store.params, scale, startCL)
      ArgsList[[i]] <- list(X, g, CnsIndx, OrdIndx, Nnorms, MaxIter, models[m],
                            store.params, scale, startCL, autoStop, ma.band, stop.tol)
      i <- i + 1
    }  # g
  }  # m
  
  ### Set up cluster
  cl <- parallel::makeCluster(cores, type = "SOCK")
  invisible(parallel::clusterEvalQ(cl, library(clustMD)))
  parallel::clusterSetRNGStream(cl)
  
  ### Run clustMD
  results <- list()
  t1 <- system.time(results <- parallel::parLapply(cl, ArgsList, clustMDlist))
  cat("Time for computation:", t1[3], "\n")
  
  ### Stop cluster
  parallel::stopCluster(cl)
  
  ### Save BIC values
  for (m in 1:length(models)) {
    for (g in 1:length(G)) {
      if (!is.null(results[[{(m - 1) * length(G) + g}]])){ 
        BICarray[g, m] <- results[[{(m - 1) * length(G) + g}]]$BIChat
        ICLarray[g, m] <- results[[{(m - 1) * length(G) + g}]]$ICLhat
      }
    }  # g
  }  # m
  
  ### Return output list and BIC matrix
  resParallel <- list(Output = results, BICarray = BICarray, ICLarray=ICLarray,
                      G = G, models = models)
  class(resParallel) <- "clustMDparallel"
  resParallel
}

#-------------------------------------------------------------------------#
