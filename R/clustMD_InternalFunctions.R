# This file contains internal clustMD functions that are not intended to
# be used explicitly by the end user.

#---------------------------------------------------------------------------#

### Function to estimate the probability that an observation is in a
### particular nominal catergory. Also estimates the first and second moment
### of the truncated normal distributiom associated with nominal variable.

#' Transforms Monte Carlo simulated data into categorical data. Calculates
#' empirical moments of latent data given categorical responses.
#' 
#' Internal function.
#'
#' @param Z a matrix of Monte Carlo simulated data.
#'
#' @return Output required for \code{clustMD} function.
#' 
#' @seealso \code{\link{clustMD}}
#' @keywords internal
#'
z.nom.diag <- function(Z) {
  # Z is a matrix of simulated vectors (each row) y is the jth column of
  # observed nominal responses
  yrep <- rep(0, dim(Z)[1])
  yrep[apply(Z, 1, max) < 0] <- 1
  yrep[yrep != 1] <- apply(Z[yrep != 1, ], 1, which.max) + 1
  
  probs <- as.vector(table(yrep)/dim(Z)[1])
  
  if (length(probs) < (dim(Z)[2] + 1)) {
    cat("ERROR:No Monte Carlo observations generating one or more levels
        of a nominal variable.", "\n", "Increasing Nnorms may solve this
        problem.")
  }
  
  Ez_nom <- matrix(NA, dim(Z)[2], dim(Z)[2] + 1)
  Ezzt_nom <- array(NA, c(dim(Z)[2], dim(Z)[2], dim(Z)[2] + 1))
  for (k in 1:(dim(Z)[2] + 1)) {
    Ez_nom[, k] <- colMeans(matrix(Z[yrep == k, ], nrow = sum(yrep == k)))
    Ezzt_nom[, , k] <- matrix(t(Z[yrep==k,])%*%Z[yrep==k, ]/sum(yrep == k),
                              dim(Z)[2], dim(Z)[2])
  }
  
  # Returns: 1. Probability of each possible response 2. Expected value of
  # latent vector for each possible response 3. Expected value of outer
  # product of latent vector for each possible response
  list(probs, Ez_nom, Ezzt_nom)
}

#---------------------------------------------------------------------------#

### Function to calculate probabilities, first and second moments of latent
### variable associated with nominal variable. Calls Z.nom.diag.

#' Calculates the first and second moments of the latent data
#' 
#' Internal function.
#'
#' @param D dimension of the latent data.
#' @param G number of mixture components.
#' @param N number of observations.
#' @param CnsIndx the number of continuous variables.
#' @param OrdIndx the sum of the number of continuous and ordinal (including
#'     binary) variables.
#' @param zlimits the truncation points for the latent data.
#' @param mu a \code{D x G} matrix of means.
#' @param Sigma a \code{D x D x G} array of covariance parameters.
#' @param Y an \code{N x J} data matrix.
#' @param J the number of variables.
#' @param K the number of levels for each variable.
#' @param norms a matrix of standard normal deviates.
#' @param nom.ind.Z the latent dimensions corresponding to each nominal
#'     variable.
#' @param patt.indx a list of length equal to the number of observed response 
#'     patterns. Each entry of the list details the observations for which that
#'     response pattern was observed.
#'
#' @return Output required for \code{clustMD} function.
#' 
#' @seealso \code{\link{clustMD}}
#' @keywords internal
#' 
z.moments <- function(D, G, N, CnsIndx, OrdIndx, zlimits, mu, Sigma, Y, J, 
                      K, norms, nom.ind.Z, patt.indx) {
  D <- J
  if (J > OrdIndx) 
    D <- OrdIndx + sum(K[(OrdIndx + 1):J] - 1)
  
  Ez.new <- array(0, c(N, D, G))
  S.new <- array(0, c(D, D, G, N))
  probs.new <- NA  # dummy required for output
  if (J > OrdIndx) {
    probs.new <- array(NA, c(J - OrdIndx, max(K[(OrdIndx + 1):J]), G))
  }
  
  # continuous
  if (CnsIndx > 0) {
      Ez.new[, 1:CnsIndx, 1:G] <- Y[, 1:CnsIndx]
      for (g in 1:G) S.new[1:CnsIndx, 1:CnsIndx, g, ] <- 
          apply(matrix(Y[, 1:CnsIndx], nrow = N), 1, vec.outer)
  }  # if
  
  for (g in 1:G) {
    # ordinal
    if (OrdIndx > CnsIndx) {
      for (p in 1:length(patt.indx)) {
        temp.moments <- dtmvnom(a = zlimits[patt.indx[[p]][1], (CnsIndx + 1):OrdIndx, 1],
                                b = zlimits[patt.indx[[p]][1], (CnsIndx + 1):OrdIndx, 2],
                                mu = mu[(CnsIndx + 1):OrdIndx, g], 
                                S = Sigma[(CnsIndx + 1):OrdIndx, (CnsIndx + 1):OrdIndx, g])
        
        Ez.new[patt.indx[[p]], (CnsIndx + 1):OrdIndx, g] <- 
          matrix(temp.moments$tmean, length(patt.indx[[p]]), OrdIndx - CnsIndx, byrow = TRUE)
        
        S.new[(CnsIndx + 1):OrdIndx, (CnsIndx + 1):OrdIndx, g, patt.indx[[p]]] <-
          (temp.moments$tvar + t(temp.moments$tvar))/2 + vec.outer(temp.moments$tmean)
      }  # p
    }  # if
    
    # Nominal
    if (J > OrdIndx) {
      for (j in (OrdIndx + 1):J) {
        Zrep <- norms[, 1:(K[j] - 1)] %*% 
          chol(Sigma[nom.ind.Z[[j - OrdIndx]], nom.ind.Z[[j - OrdIndx]], g]) +
          matrix(mu[nom.ind.Z[[j -  OrdIndx]], g], dim(norms)[1], K[j] - 1, byrow = TRUE)
        
        temp.z <- z.nom.diag(Zrep)
        probs.new[j - OrdIndx, 1:K[j], g] <- temp.z[[1]]
        
        for (k in 1:K[j]) {
        Ez.new[Y[, j] == k, nom.ind.Z[[j - OrdIndx]], g] <-
          matrix(temp.z[[2]][, k], sum(Y[, j] == k), K[j] - 1, byrow = TRUE)
        S.new[nom.ind.Z[[j - OrdIndx]], nom.ind.Z[[j - OrdIndx]], g, Y[, j] == k] <-
          array(temp.z[[3]][, , k], c(K[j] - 1, K[j] - 1, sum(Y[, j] == k)))
        }  # k
      }  # j
    }  # if
  }  # g
  list(Ez.new, probs.new, S.new)
}

#---------------------------------------------------------------------------#

### Function to calculate nominal moments and probabilities for diagonal
### models
#' Calculates the first and second moments of the latent data for diagonal models
#' 
#' Internal function.
#'
#' @param D dimension of the latent data.
#' @param G number of mixture components.
#' @param N number of observations.
#' @param CnsIndx the number of continuous variables.
#' @param OrdIndx the sum of the number of continuous and ordinal (including
#'     binary) variables.
#' @param zlimits the truncation points for the latent data.
#' @param mu a \code{D x G} matrix of means.
#' @param Sigma a \code{D x D x G} array of covariance parameters.
#' @param Y an \code{N x J} data matrix.
#' @param J the number of variables.
#' @param K the number of levels for each variable.
#' @param norms a matrix of standard normal deviates.
#' @param nom.ind.Z the latent dimensions corresponding to each nominal 
#'     variable.
#'
#' @return Output required for clustMD function.
#' @seealso \code{\link{clustMD}}
#' @keywords internal
#' 
z.moments_diag <- function(D, G, N, CnsIndx, OrdIndx, zlimits, mu, Sigma, 
  Y, J, K, norms, nom.ind.Z) {
  D <- J
  if (J > OrdIndx) 
    D <- OrdIndx + sum(K[(OrdIndx + 1):J] - 1)
  
  Ez.new <- array(NA, c(N, D, G))
  # S.new <- matrix(0, N, G)
  S2.new <- array(NA, c(D, D, G, N))
  probs.new <- NA  #dummy required for output
  if (J > OrdIndx) 
    probs.new <- array(NA, c(J - OrdIndx, max(K[(OrdIndx + 1):J]), G))
  
  for (g in 1:G) {
    # continuous
    if (CnsIndx > 0) 
      Ez.new[, 1:CnsIndx, ] <- Y[, 1:CnsIndx]
    
    # ordinal
    if (OrdIndx > CnsIndx) {
      for (i in 1:N) {
        temp.e <- truncnorm::etruncnorm(a = zlimits[i, (CnsIndx + 1):OrdIndx, 1],
                                        b = zlimits[i, (CnsIndx + 1):OrdIndx, 2],
                                        mean = mu[(CnsIndx + 1):OrdIndx, g],
                                        sd = sqrt(diag(Sigma[, , g])[(CnsIndx + 1):OrdIndx]))
        
        Ez.new[i, (CnsIndx + 1):OrdIndx, g] <- temp.e
          
        temp.v <- truncnorm::vtruncnorm(a = zlimits[i, (CnsIndx + 1):OrdIndx, 1],
                                        b = zlimits[i, (CnsIndx + 1):OrdIndx, 2],
                                        mean = mu[(CnsIndx + 1):OrdIndx, g],
                                        sd = sqrt(diag(Sigma[, , g])[(CnsIndx + 1):OrdIndx]))
        
        S2.new[(CnsIndx + 1):OrdIndx, (CnsIndx + 1):OrdIndx, g, i] <- diag(temp.v + temp.e^2,
                                                                           nrow=OrdIndx-CnsIndx) 
      }  # i
    }  # if
    
    # Nominal
    if (J > OrdIndx) {
      for (j in (OrdIndx + 1):J) {
        Zrep <- norms[, 1:(K[j] - 1)] %*% 
          chol(Sigma[nom.ind.Z[[j - OrdIndx]], nom.ind.Z[[j - OrdIndx]], g]) +
          matrix(mu[nom.ind.Z[[j - OrdIndx]], g], dim(norms)[1], K[j] - 1, byrow = TRUE)
        
        temp.z <- z.nom.diag(Zrep)
        
        probs.new[j - OrdIndx, 1:K[j], g] <- temp.z[[1]]
        
        for (k in 1:K[j]) {
          Ez.new[Y[, j] == k, nom.ind.Z[[j - OrdIndx]], g] <-
            matrix(temp.z[[2]][, k], sum(Y[, j] == k), K[j] - 1, byrow = TRUE)

          S2.new[nom.ind.Z[[j - OrdIndx]], nom.ind.Z[[j - OrdIndx]], g, Y[, j] == k] <-
            matrix(diag(diag(temp.z[[3]][, , k])), K[j] - 1, K[j] - 1, byrow = TRUE)
        }  # k
      }  # j
    }  # if
  }  # g
  list(Ez.new, probs.new, S2.new)
}

#---------------------------------------------------------------------------#

### Function to fix values for ordinal variable cutoffs based on comparing
### the sample proportions to the quantiles of the standard normal
### distribution.

#' Calculates the threshold parameters for ordinal variables.
#'
#' @param CnsIndx the number of continuous variables.
#' @param OrdIndx the sum of the number of continuous and ordinal (including
#'     binary) variables.
#' @param Y an \code{N x J} data matrix.
#' @param N number of observations.
#'
#' @return Output required for \code{clustMD} function.
#' 
#' @seealso \code{\link{clustMD}}
#' @keywords internal
#' 
perc.cutoffs <- function(CnsIndx, OrdIndx, Y, N) {
  perc.cut <- list()
  for (j in (CnsIndx + 1):OrdIndx) {
    perc.cut[[j]] <- qnorm(c(0, cumsum(table(Y[, j])/N)))
  }
  perc.cut
}

#---------------------------------------------------------------------------#

### Function to approximate the observed log likelihood
#' Approximates the observed log likelihood.
#'
#' @param N the number of observations.
#' @param CnsIndx the number of continuous variables.
#' @param G the number of mixture components.
#' @param Y an \code{N x J} data matrix.
#' @param mu a \code{D x G} matrix of means.
#' @param Sigma a \code{D x D x G} array of covariance parameters.
#' @param pi.vec the mixing weights.
#' @param patt.indx a list of length equal to the number of observed response 
#'     patterns. Each entry of the list details the observations for which that
#'     response pattern was observed.
#' @param zlimits the truncation points for the latent data.
#' @param J the number of variables.
#' @param OrdIndx the sum of the number of continuous and ordinal (including
#'     binary) variables.
#' @param probs.nom an array containing the response probabilities for each
#'     nominal variable for each cluster
#' @param model the covariance model fitted to the data.
#' @param perc.cut threshold parameters.
#' @param K the number of levels for each variable.
#'
#' @return Output required for \code{clustMD} function.
#' @import stats
#' 
#' @seealso \code{\link{clustMD}}
#' @keywords internal
#' 
ObsLogLikelihood <- function(N, CnsIndx, G, Y, mu, Sigma, pi.vec, patt.indx, 
                             zlimits, J, OrdIndx, probs.nom, model, perc.cut, K) {
  # Continuous
  logLikeCns <- rep(0, N)
  if (CnsIndx > 0) {
    densCns <- matrix(NA, N, G)
    for (g in 1:G) {
      densCns[, g] <- 
        mvtnorm::dmvnorm(matrix(Y[, 1:CnsIndx], nrow = N),
                         mean = mu[1:CnsIndx, g],
                         sigma = matrix(Sigma[1:CnsIndx, 1:CnsIndx, g], CnsIndx, CnsIndx),
                         log = TRUE)
    }
    densCns <- sweep(densCns, 2, log(pi.vec), "+")
    logLikeCns <- apply(densCns, 1, stable.probs)
  }
  
  # Categorical
  logLikeCat <- rep(0, N)
  if (J > CnsIndx) {
    # Ordinal
    densOrd <- matrix(1, N, G)
    if (OrdIndx > CnsIndx) {
      
      if(model=="BD") {
        for (g in 1:G) {
          for (p in 1:length(patt.indx)) {
            densOrd[patt.indx[[p]], g] <- 
              mvtnorm::pmvnorm(lower = zlimits[patt.indx[[p]][1], (CnsIndx + 1):OrdIndx, 1],
                               upper = zlimits[patt.indx[[p]][1], (CnsIndx + 1):OrdIndx, 2],
                               mean = mu[(CnsIndx + 1):OrdIndx, g],
                               sigma = Sigma[(CnsIndx + 1):OrdIndx, (CnsIndx + 1):OrdIndx, g])
          }  # p
        }  # g
      } else {
        OrdProbs <- array(NA, c(OrdIndx-CnsIndx, max(K[(CnsIndx+1):OrdIndx]),G))
        for(j in (CnsIndx+1):OrdIndx){
          for(g in 1:G){
            CumulProbs <- pnorm(perc.cut[[j]], mean=mu[j, g], sd=sqrt(Sigma[j, j, g]))
            for(k in 1:K[j])
              OrdProbs[j-CnsIndx, k, g] <- CumulProbs[k+1] - CumulProbs[k]
            
            densOrd[, g] <- densOrd[, g]*OrdProbs[j-CnsIndx, Y[, j], g]
          } # g
        } # j
      } # ifelse
      
    } # if
    
    # Nominal
    densNom <- matrix(1, N, G)
    if (J > OrdIndx) {
      for (g in 1:G) {
        for (j in (OrdIndx + 1):J) 
          densNom[, g] <- densNom[, g] * probs.nom[j - OrdIndx, Y[, j], g]
      }  # g
    }
    
    logDensCat <- log(densOrd) + log(densNom)
    logDensCat <- sweep(logDensCat, 2, log(pi.vec), "+")
    logLikeCat <- apply(logDensCat, 1, stable.probs)
  }
  logLike <- sum(logLikeCns + logLikeCat)
  logLike
}
#---------------------------------------------------------------------------#

### Functionto calculate the number of free parameters in each model
#' Calculates the number of free parameters for the \code{clustMD} model.
#' 
#' Internal function.
#'
#' @param model the \code{clustMD} covariance model fitted.
#' @param D the dimension of the latent data.
#' @param G the number of mixture components.
#' @param J the number of variables.
#' @param CnsIndx the number of continuous variables.
#' @param OrdIndx the sum of the number of continuous and ordinal (including
#'     binary) variables.
#' @param K a vector indicating the number of levels of each categorical 
#'     variable.
#'
#' @return Output required for \code{clustMD} function.
#' @seealso \code{\link{clustMD}}
#' @keywords internal
#' 
npars_clustMD <- function(model, D, G, J, CnsIndx, OrdIndx, K) {
  if (model == "EII") {
    npars <- 1 + (G - 1) + G * D
    # if nominal
    if (J > OrdIndx) {
      npars <- 1 + (G - 1) + G * OrdIndx + (G - 1) * (D - OrdIndx)
    }
  } else if (model == "VII") {
    npars <- G + (G - 1) + G * D
    # if nominal
    if (J > OrdIndx) {
      npars <- G + (G - 1) + (G - 1) + G * OrdIndx + (G - 1) * (D - OrdIndx)
    }
  } else if (model == "EEI") {
    npars <- 1 + (D - 1) + (G - 1) + G * D
    # if nominal
    if (J > OrdIndx) {
      npars <- 1 + (OrdIndx - 1) + (G - 1) + G * OrdIndx + (G - 1) * (D - OrdIndx)
    }
  } else if (model == "VEI") {
    npars <- G + (D - 1) + (G - 1) + G * D
    # if nominal
    if (J > OrdIndx) {
      npars <- G + (G - 1) + (OrdIndx - 1) + (G - 1) + G * OrdIndx + (G - 1) *
        (D - OrdIndx)
    }
  } else if (model == "EVI") {
    npars <- 1 + G * (D - 1) + (G - 1) + G * D
    # if nominal
    if (J > OrdIndx) {
      npars <- 1 + G * (OrdIndx - 1) + (G - 1) * (D - OrdIndx - 1) + (G - 1) +
        G * OrdIndx + (G - 1) * (D - OrdIndx)
    }
  } else if (model == "VVI") {
    npars <- G + G * (D - 1) + (G - 1) + G * D
    # if nominal
    if (J > OrdIndx) {
      npars <- G + (G - 1) + G * (OrdIndx - 1) + (G - 1) * (D - OrdIndx - 1) +
        (G - 1) + G * OrdIndx + (G - 1) * (D - OrdIndx)
    }
  } else if (model == "BD") {
  
    npars <- G * CnsIndx * (CnsIndx + 1)/2 + (G - 1) + G * D
    
    # if ordinal
    if (J > CnsIndx) {
      npars <- npars + G * {(OrdIndx - CnsIndx) * (OrdIndx - CnsIndx + 1)/2}
    }
    
    # if nominal
    if (J > OrdIndx){
      npars <- npars + G * {sum((K[(OrdIndx + 1):J] - 1) * (K[(OrdIndx + 1):J])/2)}
    }
  }
  npars
}

#---------------------------------------------------------------------------#

### Function to perform (Monte Carlo) E-step. Returns posterior
### probabilities matrix along with other precalcualted parameters for other
### functions

#' E-step of the (MC)EM algorithm
#' 
#' Internal function.
#'
#' @param N number of observations.
#' @param G number of mixture components.
#' @param D dimension of the latent data.
#' @param CnsIndx the number of continuous variables.
#' @param OrdIndx the sum of the number of continuous and ordinal (including
#'     binary) variables.
#' @param zlimits the truncation points for the latent data.
#' @param mu a D x G matrix of means.
#' @param Sigma a D x D x G array of covariance parameters.
#' @param Y an N x J data matrix.
#' @param J the number of observed variables.
#' @param K the number of levels for each variable.
#' @param norms a matrix of standard normal deviates.
#' @param nom.ind.Z the latent dimensions corresponding to each nominal 
#'     variable.
#' @param patt.indx a list of length equal to the number of observed response 
#'     patterns. Each entry of the list details the observations for which that
#'     response pattern was observed.
#' @param pi.vec mixing weights.
#' @param model the covariance model fitted to the data.
#' @param perc.cut threshold parameters.
#'
#' @return Output required for clustMD function.
#' @import stats
#' @seealso \code{\link{clustMD}}
#' @keywords internal

E.step <- function(N, G, D, CnsIndx, OrdIndx, zlimits, mu, Sigma, Y, J, K, 
                   norms, nom.ind.Z, patt.indx, pi.vec, model, perc.cut) {
  # tau
  tau.mat <- matrix(NA, N, G)
  
  # Z moments
  if (model == "BD") {
    temp.z <- z.moments(D, G, N, CnsIndx, OrdIndx, zlimits, mu, Sigma, Y, J, K,
                        norms, nom.ind.Z, patt.indx)
  } else {
    temp.z <- z.moments_diag(D, G, N, CnsIndx, OrdIndx, zlimits, mu, Sigma, Y,
                             J, K, norms, nom.ind.Z)
  }
  Ez <- temp.z[[1]]
  probs.nom <- temp.z[[2]]
  S <- temp.z[[3]]
  
  for (g in 1:G) {
    # continuous
    if (CnsIndx > 0) {
      tau.cns <- mvtnorm::dmvnorm(matrix(Y[, 1:CnsIndx], nrow = N),
                                  mean = mu[1:CnsIndx,     g], 
                                  sigma = matrix(Sigma[1:CnsIndx, 1:CnsIndx, g],
                                                 nrow = CnsIndx), log = TRUE)
    } else {
      tau.cns <- rep(0, N)
    }
    
    # ordinal
    if (OrdIndx > CnsIndx) {
      
      if (model == "BD") {
        tau.ord <- NULL
        for (p in 1:length(patt.indx)) {
          tau.ord[patt.indx[[p]]] <- 
            log(mvtnorm::pmvnorm(zlimits[patt.indx[[p]][1], (CnsIndx + 1):OrdIndx, 1],
                                 zlimits[patt.indx[[p]][1], (CnsIndx + 1):OrdIndx, 2],
                                 mean = mu[(CnsIndx + 1):OrdIndx, g],
                                 sigma = Sigma[(CnsIndx + 1):OrdIndx, (CnsIndx + 1):OrdIndx, g]))
        } # p
      } else {
        tau.ord <- matrix(NA, N, OrdIndx-CnsIndx)
        for(j in (CnsIndx+1):OrdIndx){
          resp.probs.g <- pnorm(perc.cut[[j]][2:length(perc.cut[[j]])],
                                mean=mu[j, g],
                                sd=sqrt(Sigma[j, j, g])) -
                          pnorm(perc.cut[[j]][1:(length(perc.cut[[j]])-1)],
                                mean=mu[j, g], 
                                sd=sqrt(Sigma[j, j, g]))
          
          tau.ord[, j-CnsIndx] <- resp.probs.g[Y[,j]]
        }
        tau.ord <- apply(log(tau.ord), 1, sum)
      } # ifelse
      
    } else {
      
      tau.ord <- rep(0, N)
      
    } # ifelse
    
    # nominal
    if (J > OrdIndx) {
      tau.nom <- matrix(NA, N, J - OrdIndx)
      for (j in (OrdIndx + 1):J) {
        resp.probs.g <- probs.nom[j - OrdIndx, 1:K[j], g]
        tau.nom[, j - OrdIndx] <- log(resp.probs.g[Y[, j]])
      }
      tau.nom <- rowSums(tau.nom)
    } else {
      tau.nom <- rep(0, N)
    }
    
    tau.mat[, g] <- log(pi.vec[g]) + tau.cns + tau.ord + tau.nom
  }  # g
  
  tau.new <- exp(tau.mat - apply(tau.mat, 1, stable.probs))
  
  sumTauEz.new <- apply(sweep(Ez, c(1, 3), tau.new, "*"), c(2, 3), sum)
  sumTauS.new <- apply(sweep(S, c(3, 4), t(tau.new), "*"), c(1, 2, 3), sum)
  
  # list(tau.new, sumTauEz.new, sumTauS.new, probs.nom)
  list(tau.new, sumTauEz.new, sumTauS.new, probs.nom, Ez)
}

#---------------------------------------------------------------------------#

## M-step function. Estimates mixing proportions, means and variance
## parameters

#' M-step of the (MC)EM algorithm
#' 
#' Internal function.
#'
#' @param tau a \code{N x G} matrix of cluster membership probabilities.
#' @param N number of observations.	
#' @param sumTauEz the sum across all observations of observed and expected 
#'     latent continuous values mutiplied by the posterior probability of 
#'     belonging to each cluster.
#' @param J the number of variables.
#' @param OrdIndx the sum of the number of continuous and ordinal (including
#'     binary) variables.
#' @param D dimension of the latent data.
#' @param G the number of mixture components.
#' @param Y a \code{N x J} data matrix.
#' @param CnsIndx the number of continuous variables.
#' @param sumTauS the sum across all observations of outer product of observed
#'     and expected latent continuous values mutiplied by the posterior 
#'     probability of belonging to each cluster.
#' @param model which \code{clustMD} covariance model is fitted.
#' @param a a \code{G x D} matrix of the entries of A.
#' @param nom.ind.Z the latent dimensions corresponding to each nominal variable.
#'
#' @return Output required for \code{clustMD} function.
#' 
#' @seealso \code{\link{clustMD}}
#' @keywords internal
#' 
M.step <- function(tau, N, sumTauEz, J, OrdIndx, D, G, Y, CnsIndx, sumTauS, 
                   model, a, nom.ind.Z) {
  
  if (!is.element(model, c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "BD"))) 
    print("Unknown model chosen! Choose from; EII, VII, EEI, VEI, EVI, VVI or BD")
  
  tausum <- colSums(tau)
  # Mixing weights
  pi.vec.new <- tausum/N
  
  # Mean vectors
  mu.new <- sweep(sumTauEz, 2, tausum, "/")
  # Scale mean (ONLY REQUIRED FOR NOMINALS)
  if (model != "BD") {
    if (J > OrdIndx) {
      mu.new[(OrdIndx + 1):D, ] <- 
        sweep(as.matrix(mu.new[(OrdIndx + 1):D, ]), 1, 
              rowSums(sweep(as.matrix(mu.new[(OrdIndx + 1):D, ]), 2, pi.vec.new, "*")), "-")
    }
  } 
    
  # Sigma Can probably vectorise this some way for speed
  if (model == "EII") {
    Sigma.new <- array(0, c(D, D, G))
    lambda.new <- matrix(1, G, D)
    
    # Lambda (Continuous & Ordinal)
    if (OrdIndx > 0) {
      lambda.new_A <- 
        -2 * sum(diag(t(matrix(sumTauEz[1:OrdIndx, ], nrow=OrdIndx)) %*% 
                        matrix(mu.new[1:OrdIndx, ], nrow=OrdIndx))) +
        t(diag(t(matrix(mu.new[1:OrdIndx, ], nrow=OrdIndx)) %*% mu.new[1:OrdIndx, ])) %*% tausum
      if (CnsIndx > 0) 
        lambda.new_A <- lambda.new_A + sum(diag(Y[, 1:CnsIndx] %*% t(Y[, 1:CnsIndx])))
      if (OrdIndx > CnsIndx) 
        lambda.new_A <- lambda.new_A + sum(apply(sumTauS, 3, diag)[(CnsIndx + 1):OrdIndx, ])
      
      lambda.new_A <- as.numeric(lambda.new_A)/{N * OrdIndx}
      
      for (g in 1:G) {
        Sigma.new[1:OrdIndx, 1:OrdIndx, g] <- lambda.new_A * diag(OrdIndx)
        lambda.new[g, 1:OrdIndx] <- lambda.new_A
      }
    }
    
    if (J > OrdIndx) {
      for (g in 1:G) Sigma.new[(OrdIndx + 1):D, (OrdIndx + 1):D, g] <- diag(D - OrdIndx)
    }
    
    a.new <- matrix(1, G, D)
    
  } else if (model == "VII") {
    Sigma.new <- array(0, c(D, D, G))
    lambda.new <- matrix(1, G, D)
    
    if (OrdIndx > 0) {
      lambda.new_A <- rep(NA, G)
      
      for (g in 1:G) {
        lambda.new_A[g] <- -2 * t(sumTauEz[1:OrdIndx, g]) %*% mu.new[1:OrdIndx, g] +
          t(mu.new[1:OrdIndx, g]) %*% mu.new[1:OrdIndx, g] * tausum[g]
        
        if (CnsIndx > 0) 
          lambda.new_A[g] <- lambda.new_A[g] + 
            sum(diag(tau[, g] * matrix(Y[, 1:CnsIndx], ncol=CnsIndx) %*% t(matrix(Y[, 1:CnsIndx], ncol=CnsIndx))))
        
        if (OrdIndx > CnsIndx) 
          lambda.new_A[g] <- lambda.new_A[g] + 
            sum(diag(sumTauS[(CnsIndx + 1):OrdIndx, (CnsIndx + 1):OrdIndx, g]))
        
        lambda.new_A[g] <- lambda.new_A[g]/{tausum[g] * OrdIndx}
        
        lambda.new[g, 1:OrdIndx] <- lambda.new_A[g]
        Sigma.new[1:OrdIndx, 1:OrdIndx, g] <- lambda.new_A[g] * diag(OrdIndx)
      }  # g 
    }
    
    if (J > OrdIndx) {
      lambda.new_B <- rep(NA, G)
      for (g in 1:G) {
        lambda.new_B[g] <- -2 * t(sumTauEz[(OrdIndx + 1):D, g]) %*%
          mu.new[(OrdIndx + 1):D, g] + t(mu.new[(OrdIndx + 1):D, g]) %*%
          mu.new[(OrdIndx + 1):D, g] * tausum[g] + 
          sum(diag(sumTauS[(OrdIndx + 1):D, (OrdIndx + 1):D, g]))
        
        lambda.new_B[g] <- lambda.new_B[g]/{tausum[g] * (D - OrdIndx)}
      }  # g
      
      lambda.new_B <- lambda.new_B/sum(lambda.new_B)
      
      for (g in 1:G) {
        lambda.new[g, (OrdIndx + 1):D] <- lambda.new_B[g]
        Sigma.new[(OrdIndx + 1):D, (OrdIndx + 1):D, g] <- lambda.new_B[g] * diag(D - OrdIndx)
      }  # g
    }
    
    a.new <- matrix(1, G, D)
    
  } else if (model == "EEI") {
    
    Sigma.new <- array(0, c(D, D, G))
    lambda.new <- matrix(1, G, D)
    
    # lambda (continuous & ordinal)
    if (OrdIndx > 0) {
      lambda.new_A <- -2 *
        sum(diag(t(matrix(sumTauEz[1:OrdIndx, ], nrow = OrdIndx)) %*% diag(1/a[1, 1:OrdIndx]) %*% 
                   matrix(mu.new[1:OrdIndx, ], nrow = OrdIndx))) +
        sum(diag(t(matrix(mu.new[1:OrdIndx, ], nrow = OrdIndx)) %*% diag(1/a[1, 1:OrdIndx]) %*%
                   matrix(mu.new[1:OrdIndx, ], nrow = OrdIndx)) * tausum)
      
      if (CnsIndx > 0) 
        lambda.new_A <- lambda.new_A + 
          sum(sweep(matrix(Y[, 1:CnsIndx]^2, ncol = CnsIndx), 2, a[1, 1:CnsIndx], "/"))
      if (OrdIndx > CnsIndx) 
        lambda.new_A <- lambda.new_A + 
          sum(rowSums(matrix(apply(sumTauS, 3, diag)[(CnsIndx + 1):OrdIndx, ],
                             OrdIndx - CnsIndx, G)) %*%
                (1/a[1, (CnsIndx + 1):OrdIndx]))
      
      lambda.new_A <- lambda.new_A/{N * OrdIndx}
      lambda.new[, 1:OrdIndx] <- lambda.new_A
    }
    
    # A
    a.new <- -2 * diag(sumTauEz %*% t(mu.new)) + mu.new^2 %*% tausum
    if (CnsIndx > 0) 
      a.new[1:CnsIndx] <- a.new[1:CnsIndx] +
        colSums(matrix(Y[, 1:CnsIndx]^2, N, ncol = CnsIndx))
    
    if (J > CnsIndx) 
      a.new[(CnsIndx + 1):D] <- a.new[(CnsIndx + 1):D] +
        rowSums(matrix(apply(sumTauS, 3, diag)[(CnsIndx + 1):D, ], D - CnsIndx, G))
    
    a.new <- a.new/(N * lambda.new[1, ])
    # det(A) = 1 constraint
    if (OrdIndx > 0) 
      a.new[1:OrdIndx] <- a.new[1:OrdIndx]/prod(a.new[1:OrdIndx])^(1/OrdIndx)
    
    if (J > OrdIndx) 
      a.new[(OrdIndx + 1):D] <- 1
    
    a.new <- matrix(a.new, G, D, byrow = TRUE)
    Sigma.new <- array(diag(lambda.new[1, ] * a.new[1, ], nrow = D, ncol = D), c(D, D, G))
    
  } else if (model == "VEI") {
    
    Sigma.new <- array(0, c(D, D, G))
    lambda.new <- matrix(1, G, D)
    
    # lambda (continuous & ordinal)
    if (OrdIndx > 0) {
      lambda.new_A <- -2 * 
        diag(t(matrix(sumTauEz[1:OrdIndx, ], nrow = OrdIndx)) %*% diag(1/a[1, 1:OrdIndx]) %*% matrix(mu.new[1:OrdIndx, ], nrow=OrdIndx)) +
        (t(1/a[1, 1:OrdIndx]) %*% matrix(mu.new[1:OrdIndx, ]^2, nrow = OrdIndx)) * tausum
      if (CnsIndx > 0) 
        lambda.new_A <- lambda.new_A + t(matrix(Y[, 1:CnsIndx]^2, ncol = CnsIndx) %*% 
                                           matrix(1/a[1, 1:CnsIndx], nrow = CnsIndx)) %*% tau
      if (OrdIndx > CnsIndx) 
        lambda.new_A <- lambda.new_A + (1/a[1, (CnsIndx + 1):OrdIndx]) %*% 
          matrix(apply(sumTauS, 3, diag)[(CnsIndx + 1):OrdIndx, ], nrow = OrdIndx - CnsIndx)
      lambda.new[, 1:OrdIndx] <- lambda.new_A/{OrdIndx * tausum}
    }
    
    # lambda (nominal)
    if (J > OrdIndx) {
      lambda.new_B <- -2 * 
        diag(t(matrix(sumTauEz[(OrdIndx + 1):D, ], nrow = D - OrdIndx)) %*%
               diag(1/a[1, (OrdIndx + 1):D]) %*%
               matrix(mu.new[(OrdIndx + 1):D, ], nrow=D-OrdIndx)) + 
        (t(1/a[1, (OrdIndx + 1):D]) %*% matrix(mu.new[(OrdIndx + 1):D, ]^2,
                                               nrow = D - OrdIndx)) * tausum +
        (1/a[1, (OrdIndx + 1):D]) %*% matrix(apply(sumTauS, 3, diag)[(OrdIndx + 1):D, ],
                                             nrow=D - OrdIndx)
      
      lambda.new_B <- lambda.new_B/{(D - OrdIndx) * tausum}
      lambda.new[, (OrdIndx + 1):D] <- lambda.new_B/sum(lambda.new_B)
    }
    
    # A
    a.new <- -2 * diag(sumTauEz %*% t(mu.new/t(lambda.new))) +
      t({mu.new^2/t(lambda.new)} %*% tausum)
    
    if (CnsIndx > 0) 
      a.new[1:CnsIndx] <- a.new[1:CnsIndx] + (1/lambda.new[, 1]) %*% t(tau) %*%
        matrix(Y[, 1:CnsIndx]^2, ncol = CnsIndx)
    
    if (J > CnsIndx) 
      a.new[(CnsIndx + 1):D] <- a.new[(CnsIndx + 1):D] +
        diag(matrix(apply(sumTauS, 3, diag)[(CnsIndx + 1):D, ], nrow = D - CnsIndx, ncol = G) %*%
               matrix(1/lambda.new[, (CnsIndx + 1):D], ncol = D - CnsIndx))
    
    a.new <- a.new/N
    # det(A) = 1 constraint
    if (OrdIndx > 0) 
      a.new[1:OrdIndx] <- a.new[1:OrdIndx] / prod(a.new[1:OrdIndx])^(1/OrdIndx)
    
    # identifiability constraint
    if (J > OrdIndx) 
      a.new[(OrdIndx + 1):D] <- 1
    
    a.new <- matrix(a.new, G, D, byrow = TRUE)
    Sigma.new <- array(NA, c(D, D, G))
    for (g in 1:G) Sigma.new[, , g] <-  diag(lambda.new[g, ] * a.new[1, ], nrow = D, ncol = D)
    
  } else if (model == "EVI") {
    
    Sigma.new <- array(0, c(D, D, G))
    lambda.new <- matrix(1, G, D)
    
    # lambda (continuous & ordinal)
    if (OrdIndx > 0) {
      lambda.new_A <- -2 * 
        sum(matrix(sumTauEz[1:OrdIndx, ], nrow = OrdIndx) * matrix(mu.new[1:OrdIndx, ]/t(a[, 1:OrdIndx]), OrdIndx, G)) +
        diag(t(mu.new[1:OrdIndx, ]^2) %*% (1/t(matrix(a[, 1:OrdIndx], G, OrdIndx)))) %*% tausum
      
      if (CnsIndx > 0) 
        lambda.new_A <- lambda.new_A + sum(diag(t(tau) %*% (Y[, 1:CnsIndx]^2) %*% 
          (t(matrix(1/a[, 1:CnsIndx], G, CnsIndx)))))
      
      if (OrdIndx > CnsIndx) 
        lambda.new_A <- lambda.new_A +
          sum(diag(t(matrix(apply(sumTauS, 3, diag)[(CnsIndx + 1):OrdIndx, ],
                            nrow = OrdIndx-CnsIndx)) %*%
                     t(matrix(1 / a[, (CnsIndx + 1):OrdIndx], nrow = G, ncol = OrdIndx - CnsIndx))))
      
      lambda.new_A <- lambda.new_A / {OrdIndx * N}
      lambda.new[, 1:OrdIndx] <- lambda.new_A
    }
    
    # A
    a.new <- t(-2 * sumTauEz * mu.new + sweep(mu.new^2, 2, tausum, "*"))
    if (CnsIndx > 0) 
      a.new[, 1:CnsIndx] <- a.new[, 1:CnsIndx] + t(tau) %*% (Y[, 1:CnsIndx]^2)
    
    if (J > CnsIndx) 
      a.new[, (CnsIndx + 1):D] <- a.new[, (CnsIndx + 1):D] +
        t(apply(sumTauS, 3, diag)[(CnsIndx + 1):D, ])
    
    a.new <- a.new/sweep(lambda.new, 1, tausum, "*")
    
    # det(A) = 1 constraint
    if (OrdIndx > 0) 
      a.new[, 1:OrdIndx] <- 
        sweep(matrix(a.new[, 1:OrdIndx], G, OrdIndx), 1, 
              apply(matrix(a.new[, 1:OrdIndx], G, OrdIndx), 1, prod)^(1/OrdIndx), "/")
    
    # scale
    if (J > OrdIndx) 
      a.new[, (OrdIndx + 1):D] <- 
        sweep(matrix(a.new[, (OrdIndx + 1):D], G, D - OrdIndx), 2, 
              colSums(matrix(a.new[, (OrdIndx + 1):D], G, D - OrdIndx)), "/")
    
    Sigma.new <- array(NA, c(D, D, G))
    
    for (g in 1:G) Sigma.new[, , g] <- lambda.new[g, ] * diag(a.new[g, ])
    
  } else if (model == "VVI") {
    
    Sigma.new <- array(0, c(D, D, G))
    lambda.new <- matrix(1, G, D)
    
    # lambda (continuous & ordinal)
    if (OrdIndx > 0) {
      lambda.new_A <- -2 * 
        diag(t(matrix(sumTauEz[1:OrdIndx, ], nrow=OrdIndx)) %*% 
               matrix(mu.new[1:OrdIndx, ]/t(a[, 1:OrdIndx]), OrdIndx, G)) +
        diag(t(matrix(mu.new[1:OrdIndx, ]^2, nrow=OrdIndx)) %*%
               (1/t(matrix(a[, 1:OrdIndx], G, OrdIndx)))) * tausum
      
      if (CnsIndx > 0) 
        lambda.new_A <- lambda.new_A + diag(t(tau) %*% matrix(Y[, 1:CnsIndx]^2, ncol = CnsIndx) %*%
                                              (t(matrix(1/a[, 1:CnsIndx], G, CnsIndx))))
      if (OrdIndx > CnsIndx) 
        lambda.new_A <- lambda.new_A + 
          diag(t(matrix(apply(sumTauS, 3, diag)[(CnsIndx + 1):OrdIndx, ],
                        nrow = OrdIndx - CnsIndx)) %*% t(matrix(1/a[, (CnsIndx + 1):OrdIndx],
                                                                nrow = G, ncol = OrdIndx - CnsIndx)))
      
      lambda.new_A <- lambda.new_A / {OrdIndx * tausum}
      lambda.new[, 1:OrdIndx] <- lambda.new_A
    }
    
    # lambda (nominal)
    if (J > OrdIndx) {
      lambda.new_B <- -2 * 
        diag(t(sumTauEz[(OrdIndx + 1):D, ]) %*% 
               matrix(mu.new[(OrdIndx + 1):D, ]/t(a[, (OrdIndx + 1):D]), D - OrdIndx, G)) +
        diag(t(mu.new[(OrdIndx + 1):D, ]^2) %*% (1/t(matrix(a[, (OrdIndx + 1):D],
                                                            G, D - OrdIndx)))) *
        tausum + diag(t(apply(sumTauS, 3, diag)[(OrdIndx + 1):D, ]) %*%
                        t(matrix(1/a[, (OrdIndx + 1):D], G, D - OrdIndx)))
      
      lambda.new_B <- lambda.new_B/{(D - OrdIndx) * tausum}
      lambda.new[, (OrdIndx + 1):D] <- lambda.new_B/sum(lambda.new_B)
    }
    
    # A
    a.new <- t(-2 * sumTauEz * mu.new + sweep(mu.new^2, 2, tausum, "*"))
    
    if (CnsIndx > 0) 
      a.new[, 1:CnsIndx] <- a.new[, 1:CnsIndx] + t(tau) %*% matrix(Y[, 1:CnsIndx]^2, ncol = CnsIndx)
    
    if (J > CnsIndx) 
      a.new[, (CnsIndx + 1):D] <- a.new[, (CnsIndx + 1):D] +
        t(apply(sumTauS, 3, diag)[(CnsIndx + 1):D, ])
    
    a.new <- a.new/sweep(lambda.new, 1, tausum, "*")
    # det(A) = 1 constraint
    if (OrdIndx > 0) 
      a.new[, 1:OrdIndx] <- 
        sweep(matrix(a.new[, 1:OrdIndx], G, OrdIndx), 1, 
              apply(matrix(a.new[, 1:OrdIndx], G, OrdIndx), 1, prod)^(1/OrdIndx), "/")
    
    # scale
    if (J > OrdIndx) 
      a.new[, (OrdIndx + 1):D] <- 
        sweep(matrix(a.new[, (OrdIndx + 1):D], G, D - OrdIndx), 2, 
              colSums(matrix(a.new[, (OrdIndx + 1):D], G, D - OrdIndx)), "/")
    
    for (g in 1:G) 
      Sigma.new[, , g] <- lambda.new[g, ] * diag(a.new[g, ])
    
  } else if (model == "BD") {
    
    Sigma.new <- array(0, c(D, D, G))
    B <- array(0, c(D, D, G))
    a.new <- matrix(NA, G, D)
    lambda.new <- matrix(NA, G, D)
    
    for (g in 1:G) {
      B[, , g] <- sumTauS[, , g] - 2 * sumTauEz[, g] %*% t(mu.new[, g]) + tausum[g] *
        vec.outer(mu.new[, g])
      
      # Continuous
      if (CnsIndx > 0) 
        {
          Sigma.new[1:CnsIndx, 1:CnsIndx, g] <- B[1:CnsIndx, 1:CnsIndx, g] / tausum[g]
        }  # if
      
      # Ordinal
      if (OrdIndx > CnsIndx) {
        Sigma.new[(CnsIndx + 1):OrdIndx, (CnsIndx + 1):OrdIndx, g] <- 
          B[(CnsIndx + 1):OrdIndx, (CnsIndx + 1):OrdIndx, g] / tausum[g]
      }  # if
      
      # Nominal
      if (J > OrdIndx) {
        for (j in 1:(J - OrdIndx)) {
          Sigma.new[nom.ind.Z[[j]], nom.ind.Z[[j]], g] <- 
            B[nom.ind.Z[[j]], nom.ind.Z[[j]], g] / tausum[g]
        }  # j
      }  # if
      
    }  # g
    
  } else {
    print("Unknown model chosen! Choose from; EII, VII, EEI, VEI, EVI, VVI or BD")
  }
  
  list(pi.vec.new, mu.new, lambda.new, a.new, Sigma.new)
}

#---------------------------------------------------------------------------#
