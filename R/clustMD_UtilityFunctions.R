# This file contains some utility functions that I have written. These may
# be useful outside of the clustMD package.

#-------------------------------------------------------------------------------#

### Function returns outer product of a vector with itself
#' Calculate the outer product of a vector with itself
#' 
#' This function calculates the outer product of a vector with itself.
#'
#' @param x a numeric vector.
#'
#' @return Returns the outer product of the vector \code{x} with itself.
#'
#' @keywords internal
#' 
vec.outer <- function(x) {
  # applying this function to a matrix produces a matrix where each column
  # is the stacked columns of the outer product matrix
  x %o% x
}

#-------------------------------------------------------------------------------#

### Function to check if response pattern is the same.  returns TRUE/FALSE
### depending on whether patterns are equal or not.
#' Check if response patterns are equal
#' 
#' Checks whether response patterns are equal or not and returns \code{TRUE}
#' or \code{FALSE} reprectively.
#'
#' @param x a numeric vector.
#' @param patt a vector to compare \code{x} to.
#'
#' @return Returns \code{TRUE} if \code{x} and \code{patt} are exactly the same
#'     and \code{FALSE} otherwise.
#' 
#' @note Used internally in clustMD function.
#' @keywords internal
#' 

patt.equal <- function(x, patt) {
  if (length(x) != length(patt)) 
    print("patt.equal: Vectors have different lengths")
  
  temp <- sum(x == patt)
  if (temp == length(x)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#-------------------------------------------------------------------------------#

### Function to perform stable calculation of log denominator for tau
### Computes the log of a sum when supplied with individual logs
#' Stable computation of the log of a sum
#' 
#' Function takes a numeric vector and returns the log of the sum of the
#' elements of that vector. Calculations are done on the log scale for 
#' stability.
#'
#' @param s a numeric vector.
#'
#' @return The log of the sum of the elements of \code{s}
#' @keywords internal
#'
stable.probs <- function(s) {
  s.max <- max(s)
  indx <- which.max(s)
  alpha <- s.max + log(1 + sum(exp(s[-indx] - s.max)))
  alpha
}

#-------------------------------------------------------------------------------#

### Parallel coordinates plot for cluster means adapted for clustMD
#' Parallel coordinates plot adapted for \code{clustMD} output
#' 
#' Produces a parallel coordinates plot as \code{parcoord} in the \code{MASS}
#' library with some minor adjustments.
#'
#' @param x a matrix or data frame who columns represent variables. Missing 
#'     values are allowed.
#' @param col a vector of colours, recycled as necessary for each observation.
#' @param xlabels a character vector of variable names for the x axis.
#' @param lty a vector of line types, recycled as necessary for each
#'     observation.
#' @param var.label if TRUE, each variable's axis is labelled with maximum and
#'     minimum values.
#' @param xlab label for the X axis.
#' @param ylab label for the Y axis.
#' @param ... further graphics parameters which are passed to \code{matplot}.
#'
#' @return A parallel coordinates plot is drawn with one line for each cluster.
#' @references Wegman, E. J. (1990) Hyperdimensional data analysis using 
#'     parallel coordinates. Journal of the American Statistical Association 
#'     85, 664-675.
#' 
#'     Venables, W. N. and Ripley, B. D. (2002) Modern Applied Statistics with
#'     S. Fourth edition. Springer.
#' 
#' @import graphics
#' 
#' @keywords device internal
#' 

clustMDparcoord <- function(x, col = 1, xlabels=NULL, lty = 1, var.label = FALSE,
                            xlab = "", ylab = "", ...) {
  par(mar = c(4.5, 3, 1.5, 4), xpd = TRUE)
  if(is.null(xlabels)) xlabels <- substr(colnames(x), 1, 7)
  rx <- apply(x, 2L, range, na.rm = TRUE)
  x <- apply(x, 2L, function(x) (x - min(x, na.rm = TRUE))/(max(x, na.rm = TRUE) - 
    min(x, na.rm = TRUE)))
  matplot(1L:ncol(x), t(x), type = "l", col = col, lty = lty, xlab = xlab, 
    ylab = ylab, axes = FALSE, lwd = 2.5, ...)
  legend("topright", inset = c(-0.13, 0), legend = 1:nrow(x), fill = col, 
    col = col, title = "Cluster", bty = "n")
  axis(1, at = 1L:ncol(x), labels = xlabels, line = 0.5, 
    las = 2, cex.axis=0.8)
  for (i in 1L:ncol(x)) {
    lines(c(i, i), c(0, 1), col = "grey70")
    if (var.label) 
      text(c(i, i), c(0, 1), labels = format(rx[, i], digits = 3), xpd = NA, 
        offset = 0.3, pos = c(1, 3), cex = 0.7)
  }
  invisible()
}

#-------------------------------------------------------------------------------#

# % % dtmvnmom.R Date: 10/6/2016 [Kan & Robotti (2016)] % This program
# computes the mean and covariance matrix of Y, % where Y is a doubly
# truncated multivariate normal with parameters % mu and S, and lower and
# upper truncation limits of a and b.  % Input: % a: lower truncation
# limits % b: upper truncation limits % mu: an nx1 vector of mean of X %
# S: an nxn covariance matrix of X % Output % muY: E[Y] = E[X|a<X<b] %
# varY: Var[Y] = Var[X|a<X<b] %

#-------------------------------------------------------------------------------#


#' Helper internal function for \code{dtmvnom()}
#' 
#' Internal function.
#'
#' @param a a vector of lower thresholds.
#' @param b a vector of upper thresholds.
#' @param S the covariance matrix of the untruncated distribution.
#'
#' @return Output required for dtmvnom function.
#' @references Kan, R., & Robotti, C. (2016). On Moments of Folded and 
#'     Truncated Multivariate  Normal Distributions. Available at SSRN.
#' @import stats
#' @keywords distribution internal
#' 
qfun <- function(a, b, S) {
  n = length(a)
  s = sqrt(diag(S))
  
  if (n == 1) {
    qa = dnorm(a/s)/s
    qb = dnorm(b/s)/s
    return(list(qa = qa, qb = qb))
  } else {
    qa = matrix(0, n, 1)
    qb = matrix(0, n, 1)
    for (i in 1:n) {
      ind = 1:n
      ind = ind[-i]
      B = matrix(S[ind, i], ncol = 1)/S[i, i]
      RR = S[ind, ind] - B %*% matrix(S[i, ind], nrow = 1)
      a1 = a[ind]
      b1 = b[ind]
      
      if (is.finite(a[i])) 
        {
          # if(n==5){ qa[i] = dnorm(a[i], 0, s[i])*quadncdfab(a1-B*a[i], b1-B*a[i],
          # RR) }else{
          qa[i] = dnorm(a[i], 0, s[i]) * 
            mvtnorm::pmvnorm(a1, b1, as.vector(B * a[i]), sigma = RR)
          # } # ifelse
        }  # if
      
      if (is.finite(b[i])) 
        {
          # if(n==5){ qb[i] = dnorm(b[i], 0, s[i])*quadncdfab(a1-B*b[i], b1-B*b[i],
          # RR) }else{
          qb[i] = dnorm(b[i], 0, s[i]) * 
            mvtnorm::pmvnorm(a1, b1, as.vector(B * b[i]), sigma = RR)
          # } # ifelse
        }  # if
    }  # i
    return(list(qa = qa, qb = qb))
  }  # ifelse
}

#-------------------------------------------------------------------------------#


#' Return the mean and covariance matrix of a truncated multivariate normal
#' distribution
#' 
#' This function returns the mean and covariance matrix of a truncated 
#' multivariate normal distribution. It takes as inputs a vector of lower 
#' thresholds and another of upper thresholds along with the mean and 
#' covariance matrix of the untruncated distribution. This function follows
#' the method proposed by Kan \& Robotti (2016).
#'
#' @param a a vector of lower thresholds.
#' @param b a vector of upper thresholds.
#' @param mu the mean of the untruncated distribution.
#' @param S the covariance matrix of the untruncated distribution.
#' 
#' @return Returns a list of two elements. The first element, \code{tmean}, is
#'     the mean of the truncated multivariate normal distribution. The second 
#'     element, \code{tvar}, is the covariance matrix of the truncated 
#'     distribution.
#'     
#' @references Kan, R., & Robotti, C. (2016). On Moments of Folded and 
#'     Truncated Multivariate  Normal Distributions. Available at SSRN.
#' 
#' @import stats
#' @keywords distribution internal

dtmvnom <- function(a, b, mu, S) {
  # transposing for some reason?? May not be necessary in R if(length(a)>1)
  # a = t(a) if(length(b)>1) b = t(b) if(length(mu)>1) mu = t(mu)
  
  # S must be a matrix
  S <- as.matrix(S)
  
  n = length(mu)
  s = sqrt(diag(S))
  
  if (n == 1) {
    
    a1 = (a - mu)/s
    b1 = (b - mu)/s
    p = pnorm(b1) - pnorm(a1)
    muY = mu + (dnorm(a1) - dnorm(b1)) * s / p
    
    if (is.infinite(a)) {
      a = 0
    }
    
    if (is.infinite(b)) {
      b = 0
    }
    
    varY = S + (mu - muY) * muY + (a * dnorm(a1) - b * dnorm(b1)) * s / p
    return(list(tmean = as.vector(muY), tvar = varY))
    
  } else {
    
    a1 = a - mu
    b1 = b - mu
    
    # if(n==4){ p = quadncdfab(a1,b1,S) }else{
    p = mvtnorm::pmvnorm(a, b, mu, sigma = S)
    # } # if
    
    # check this line
    temp1 = qfun(a1, b1, S)
    qa <- temp1[[1]]
    qb <- temp1[[2]]
    q = qa - qb
    
    muY = mu + S %*% q/p
    D = matrix(0, n, n)
    
    check_a <- is.finite(a)
    check_b <- is.finite(b)
    for (i in 1:n) {
      if (check_a[i]) 
        D[i, i] = a[i] * qa[i]
      
      if (check_b[i]) 
        D[i, i] = D[i, i] - b[i] * qb[i]
      
      # ?? go over from here down
      ind = 1:n
      ind <- ind[-i]
      RR = S[ind, ind] - matrix(S[ind, i], ncol = 1) %*% 
        matrix(S[i, ind], nrow = 1)/S[i, i]
      
      if (is.infinite(a[i])) {
        wa = matrix(0, n - 1, 1)
      } else {
        ma = mu[ind] + (S[ind, i]/S[i, i]) * a1[i]
        # ma=0
        temp2 = qfun(a[ind] - ma, b[ind] - ma, RR)
        qa1 <- temp2[[1]]
        qb1 <- temp2[[2]]
        wa = qa[i] * ma + dnorm(a[i], mu[i], s[i]) * RR %*% (qa1 - qb1)
      }
      
      if (is.infinite(b[i])) {
        wb = matrix(0, n - 1, 1)
      } else {
        mb = mu[ind] + (S[ind, i]/S[i, i]) * b1[i]
        temp2 = qfun(a[ind] - mb, b[ind] - mb, RR)
        qa2 <- temp2[[1]]
        qb2 <- temp2[[2]]
        wb = qb[i] * mb + dnorm(b[i], mu[i], s[i]) * RR %*% (qa2 - qb2)
      }
      D[i, ind] = wa - wb
    }  # i
    varY = S + S %*% (D - q %*% t(muY))/p
    return(list(tmean = as.vector(muY), tvar = varY))
  }  # ifelse
}

#------------------------------------------------------------------------------#

# modal function
#' Calculate the mode of a sample
#'
#' @param x a vector containing the sample values.
#'
#' @return The mode of the sample. In the case of a tie, the minimum is 
#'     returned.
#'
#' @keywords internal
modal.value = function(x)
{
  s = 0
  r = 0
  mde = 0
  x = sort(x)
  y = unique(x)
  for(i in 1:length(y))
  {
    s = length(which(x==y[i]))
    if(s>r){
      mde = y[i]
      r = s
    }else if(s==r){
      mde = c(mde,y[i])
    }
  }
  min(mde)
}

### Function to extract relevant output from the list returned by
### clustMDparallel resParallel is the entire object returned by
### clustMDparallel nClus and covModel are the number of clusters and the
### covariance model desired

#' Extracts relevant output from \code{clustMDparallel} object
#' 
#' This function takes a \code{clustMDparallel} object, a number of clusters
#' and a covariance model as inputs. It then returns the output corresponding
#' to that model. If the particular model is not contained in the 
#' \code{clustMDparallel} object then the function returns an error.
#'
#' @param resParallel a \code{clustMDparallel} object.
#' @param nClus the number of clusters in the desired output.
#' @param covModel the covariance model of the desired output.
#'
#' @return A \code{clustMD} object containing the output for the relevant 
#'     model.
#' @export
#'
getOutput_clustMDparallel <- function(resParallel, nClus, covModel) {
  # checks
  if (!is.element(covModel, resParallel$models)) {
    stop("covModel is not among the models fitted!")
  }
  if (!is.element(nClus, resParallel$G)) {
    stop("nClus is not among the models fitted!")
  }
  
  ### Cluster and model indices
  g <- which(resParallel$G == nClus)
  m <- which(resParallel$models == covModel)
  
  modelOut <- resParallel$Output[[((m - 1) * length(resParallel$G) + g)]]
}

#------------------------------------------------------------------------------#