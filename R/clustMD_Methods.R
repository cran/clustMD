# clustMD and clustMDparallel return objects of class 'clustMD' and
# 'clustMDparallel' respectively. This file contains some basic methods
# for these classes of objects

#######################################################
##### Methods for single clustMD output named res #####
#######################################################

##### PLOTTING METHOD FOR CLUSTMD OBJECTS This function produces plots for
##### clustMD objects. The first plot is a parallel coordinates plot of the
##### cluster means, the second batch of plots are either barplots or heatmaps
##### of covariance matrices (depending on the covariance model) and lastly a
##### histogram of clustering uncertainties is produced.

#' Plotting method for objects of class \code{clustMD}
#' 
#' Plots a parallel coordinates plot and dot plot of the estimated cluster 
#' means, a barplot of the variances by cluster for diagonal covariance models
#' or a heatmap of the covariance matrix for non-diagonal covariance 
#' structures, and a histogram of the clustering uncertainties for each 
#' observation.
#'
#' @param x a \code{clustMD} object.
#' @param ... further arguments passed to or from other methods.
#'
#' @return Prints graphical summaries of the fitted model as detailed above.
#' @export
#'
#' @references McParland, D. and Gormley, I.C. (2016). Model based clustering 
#'     for mixed data: clustMD. Advances in Data Analysis and Classification, 
#'     10 (2):155-169.
#'     
#' @import ggplot2
#' @import graphics
#' @seealso \code{\link{clustMD}}
#' @keywords device
#' 
plot.clustMD <- function(x, ...) {
  .pardefault <- par(no.readonly=TRUE)

  ## Distinct colour palette
  clustMDcols <- c("gray26", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
                   "#0072B2", "#D55E00", "#CC79A7", "red3", "green3",
                   "slateblue", "darkorange", "skyblue1", "violetred4",
                   "forestgreen", "steelblue4", "slategrey", "brown",
                   "darkseagreen", "olivedrab",  "royalblue", "tomato4",
                   "cyan2", "springgreen2")
  
  grDevices::palette(clustMDcols)
  cols_G <- 1:x$G
    
  par(ask = TRUE)
  ## Convergence assessment
  if(!is.null(x$paramlist)) {
    likeStore <- na.omit(x$likelihood.store)
    likeData <- as.data.frame(likeStore)
    
    likePlotTitle <- "EM Convergence"
    if(ncol(x$Y) < ncol(x$Sigma)) likePlotTitle <- "MCEM Convergence"
    likePlot <- ggplot(data=likeData, aes(x <- 1:length(likeStore), y=likeStore)) +
      theme_bw()
    likePlot <- likePlot + geom_point(col="darkcyan") + geom_line(col="darkcyan")
    likePlot <- likePlot +
      labs(title=likePlotTitle, x="Iteration", y="Estimated Log Likelihood")
    print(likePlot)
  }
  
  ## Parallel coordinates plot of cluster means
  clustMDparcoord(t(x$means), xlabels <- x$VarNames_sht, col = cols_G,
                  var.label = T, ylab = expression(hat(mu)[g]), cex.lab = 1.3,
                  mgp = c(0.5, 0, 0))
    
  ## Dot plot of cluster means
  Mu <- x$means
  rownames(Mu) <- x$VarNames_sht
  Mu <- reshape2::melt(Mu)
  Mu$Var2 <- as.factor(Mu$Var2)
  mean_dot <- ggplot(Mu, aes_string(x = 'Var1', y='value')) +
    geom_point(col=clustMDcols[cols_G[Mu$Var2]], shape=as.character(Mu$Var2), size=5)
  mean_dot <- mean_dot + coord_flip() +
    labs(y=expression(hat(mu)[g]), x=element_blank(), title="Estimated Cluster Means")
  mean_dot <- mean_dot + theme_bw()
  print(mean_dot)
  
  ## Barchart/heatmap of covariance matrices depending on model
  if (x$model != "BD") {
    
    # for diagonal models, extract the diagonals of the covariance matrices
    VarMat <- apply(x$Sigma, 3, diag)
    rownames(VarMat) <- x$VarNames_sht
    VarMat <- 
      reshape2::melt(VarMat, varnames=c("Variable", "Cluster"), value.name="Variance")
    VarMat$Cluster <- as.factor(VarMat$Cluster)
    
    VarBar <- ggplot(data=VarMat, aes_string(x='Variable', y='Variance', fill='Cluster')) +
      geom_bar(stat="identity", position=position_dodge(), colour="black") +
      scale_fill_manual(values=clustMDcols[1:x$G])
    VarBar <- VarBar + theme_bw() + labs(title="Estimated Cluster Variances")
    VarBar <- VarBar + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5), 
            axis.title.x = element_blank())
    print(VarBar)
      
  } else {
    
    tempSigma <- x$Sigma
    colnames(tempSigma) <- x$VarNames_sht
    rownames(tempSigma) <- x$VarNames_sht
    # reorder columns so plot is the correct way around
    tempSigma <- tempSigma[, nrow(tempSigma):1 ,]
    
    longData <- reshape2::melt(tempSigma[, , 1], value.name = "Covariance")
    covPlot <- ggplot(longData,
                  aes_string(x = 'Var1', y = 'Var2', fill = 'Covariance'))
    covPlot <- covPlot + geom_tile()
    covPlot <- covPlot + viridis::scale_fill_viridis()
    covPlot <- covPlot + scale_x_discrete(expand = c(0, 0), position = "top")
    covPlot <- covPlot + scale_y_discrete(expand = c(0, 0))
    covPlot <- covPlot + coord_equal()
    covPlot <- covPlot + theme_bw() + theme(legend.title=element_blank())
    covPlot <- covPlot + ggtitle("Covariance Matrix: Cluster 1")
    covPlot <- covPlot + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1), 
            axis.title.x = element_blank(), axis.title.y = element_blank())
    print(covPlot)
    
    if(x$G > 1){
      for(g in 2:x$G){
        longData <- reshape2::melt(tempSigma[, , g], value.name = "Covariance")
        
        covPlot <- ggplot(longData,
                      aes_string(x = 'Var1', y = 'Var2', fill = 'Covariance'))
        covPlot <- covPlot + geom_tile()
        covPlot <- covPlot + viridis::scale_fill_viridis()
        covPlot <- covPlot + scale_x_discrete(expand = c(0, 0), position = "top")
        covPlot <- covPlot + scale_y_discrete(expand = c(0, 0))
        covPlot <- covPlot + coord_equal()
        covPlot <- covPlot + theme_bw() + theme(legend.title=element_blank())
        covPlot <- covPlot + ggtitle(paste("Covariance Matrix: Cluster ", g, sep=""))
        covPlot <- covPlot + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1), 
                axis.title.x = element_blank(), axis.title.y = element_blank())
        print(covPlot)
      }
    }
  }  # ifelse
    
  ## Clustering uncertainty
  uncertainty <- 1 - apply(x$tau, 1, max)
  uncertain_dat = as.data.frame(uncertainty)
  uHist <- ggplot(data=uncertain_dat, aes(uncertainty))  + theme_bw()
  uHist <- uHist + 
    geom_histogram(breaks=seq(0, 1-(1/x$G), by = 0.025), col="black", fill="darkcyan")
  uHist <- uHist + 
    labs(title="Histogram of Clustering Uncertainty", 
         x = paste("Clustering Uncertainty \n (max. possible uncertainty for ",
                   x$G, " clusters is ", round(1-(1/x$G), 3), ")", sep=""),
         y = "Frequency")
  print(uHist)
  
  par(.pardefault)
  grDevices::palette("default")
}  # plot.clustMD

#------------------------------------------------------------------------------#

##### BASIC PRINT FUNCTION FOR CLUSTMD OBJECTS This basic function prints to
##### screen the number of clusters and the model fitted to the data. It also
##### prints the associated BIC value.
#' Print basic details of \code{clustMD} object.
#' 
#' Prints a short summary of a \code{clustMD} object to screen. Details the
#' number of clusters fitted as well as the covariance model and the estimated
#' BIC.
#'
#' @param x a \code{clustMD} object.
#' @param ... further arguments passed to or from other methods.
#'
#' @return Prints summary details, as described above, to screen.
#' @export
#'
#' @seealso \code{\link{clustMD}}
#' 
#' @keywords print 
#' 
print.clustMD <- function(x, ...) {
  cat("The number of mixture components is", x$G, "and the covariance model is",
      x$model, "\nThe estimated BIC for this model is", x$BIChat)
}

#------------------------------------------------------------------------------#

##### BASIC SUMMARY FUNCTION FOR CLUSTMD OBJECTS This basic summary function
##### produces the output from print above along with a table of the cluster
##### membership and a matrix detailing the cluster means.
#' Summarise \code{clustMD} object
#' 
#' Prints a summary of a \code{clustMD} object to screen. Details the number
#' of clusters fitted as well as the covariance model and the estimated BIC.
#' Also prints a table detailing the number of observations in each cluster and
#' a matrix of the cluster means.
#'
#' @param object a \code{clustMD} object.
#' @param ... further arguments passed to or from other methods.
#'
#' @return Prints summary of \code{clustMD} object to screen, as detailed above.
#' @export
#'
#' @seealso \code{\link{clustMD}}
#' @keywords print
#' 
summary.clustMD <- function(object, ...) {
  # basic fit details
  print.clustMD(object)
  cat("\n")
    
  # table of cluster membership
  cat("\nThe number of observations in each cluster is:")
  print(table(object$cl))
    
  # cluster means
  cat("\nThe estimated means of the observed and latent dimensions by cluster:\n")
  print(object$means)
}

#------------------------------------------------------------------------------#

###############################################################
##### Methods for clustMDparallel object named x ####
###############################################################

##### BASIC PRINT FUNCTION FOR CLUSTMDPARALLEL OBJECTS This basic function
##### prints to screen the numbers of clusters and the models fitted to the
##### data. It also states which model was best according to the estimated
##### BIC.

#' Print basic details of \code{clustMDparallel} object
#' 
#' Prints basic details of \code{clustMDparallel} object. Outputs the different
#' numbers of clusters and the different covariance structures fitted to the
#' data. It also states which model was optimal according to the estimated BIC
#' criterion.
#'
#' @param x a \code{clustMDparallel} object.
#' @param ... further arguments passed to or from other methods.
#'
#' @return Prints details described above to screen.
#' @export
#'
#' @seealso \code{\link{clustMD}}
#' @keywords print
print.clustMDparallel <- function(x, ...) {
  bestMod <- which(x$BICarray == max(x$BICarray), arr.ind = TRUE)
  cat("The numbers of mixture components fitted were:\n", x$G, 
      "\nThe covariance models fitted were:\n", x$models)
  cat("\n")
  cat("\nThe best fitting model had", x$G[bestMod[1]], "components and a", 
      x$models[bestMod[2]], "covariance structure.")
}

#------------------------------------------------------------------------------#

##### BASIC SUMMARY FUNCTION FOR CLUSTMDPARALLEL OBJECTS This basic summary
##### function produces the output from print above along with a table of the
##### cluster membership and a matrix detailing the cluster means for the best
##### fitting model.

#' Prints a summary of a clustMDparallel object to screen.
#' 
#' Prints the different numbers of clusters and covariance models fitted and
#' indicates the optimal model according to the estimated BIC criterion. The
#' estimated BIC for the optimal model is printed to screen along with a table
#' of the cluster membership and the matrix of cluster means for this optimal
#' model.
#'
#' @param object a \code{clustMDparallel} object.
#' @param ... further arguments passed to or from other methods.
#'
#' @return Prints a summary of the \code{clustMDparallel} object to screen, as
#'     detailed above.
#' @export
#'
#' @seealso \code{\link{clustMD}}
#' @keywords print
#' 
summary.clustMDparallel <- function(object, ...) {
  # basic fit details
  print.clustMDparallel(object)
  cat("\n")
    
  bestMod <- which(object$BICarray == max(object$BICarray), arr.ind = TRUE)
  bestRes <- getOutput_clustMDparallel(object, nClus = object$G[bestMod[1]], 
        covModel = object$models[bestMod[2]])
  cat("\nBest model:\n")
  summary.clustMD(bestRes)
}

#------------------------------------------------------------------------------#

##### PLOTTING METHOD FOR CLUSTMDPARALLEL OBJECTS This function produces plots
##### for clustMDparallel objects. The first plot is a line plot of the
##### estimated BIC values for the models fitted. The plots for the best
##### fitting model are then produced as in plot.clustMD()

#' Summary plots for a clustMDparallel object
#' 
#' Produces a line plot of the estimated BIC values corresponding to each 
#' covariance model against the number of clusters fitted. For the optimal model
#' according to this criteria, a parallel coordinates plot of the cluster means
#' is produced along with a barchart or heatmap of the covariance matrices for
#' each cluster and a histogram of the clustering uncertainties.
#'
#' @param x a \code{clustMDparallel} object.
#' @param ... further arguments passed to or from other methods.
#'
#' @return Produces a number of plots as detailed above.
#' @import ggplot2
#' @import graphics
#' @export
#'
plot.clustMDparallel <- function(x, ...) {
  .pardefault <- par(no.readonly=TRUE)
  
  clustMDcols <- c("gray26", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
                   "#0072B2", "#D55E00", "#CC79A7", "red3", "green3",
                   "slateblue", "darkorange", "skyblue1", "violetred4",
                   "forestgreen", "steelblue4", "slategrey", "brown",
                   "darkseagreen", "olivedrab",  "royalblue", "tomato4",
                   "cyan2", "springgreen2")
  
  grDevices::palette(clustMDcols)
  
  cols_mod <- 1:length(x$models)
  GGcols_mod <- clustMDcols[rep(cols_mod, rep(max(x$G), length(x$models)))]  

  ## BIC plot
  BICdat <- reshape2::melt(x$BICarray) 
  BICdat$Var1 <- as.numeric(BICdat$Var1)
  BICplot <- ggplot(data=BICdat, aes_string(x='Var1', y='value', group='Var2',
                                            colour='Var2'))
  BICplot <- BICplot + geom_line() + geom_point() +
    scale_colour_manual(values=clustMDcols[cols_mod])
  BICplot <- BICplot +
    labs(y=expression(widehat("BIC")), x="G") + ggtitle("Model Selection")
  BICplot <- BICplot + theme_bw() + theme(legend.title=element_blank()) 
  print(BICplot)
  
  ## Get output for the optimal model
  bestMod <- which(x$BICarray == max(x$BICarray), arr.ind = TRUE)
  bestRes <- getOutput_clustMDparallel(x, nClus = x$G[bestMod[1]], 
                                       covModel = x$models[bestMod[2]])
  
  ## Remaining plots for best models given by plot.clustMD
  plot.clustMD(bestRes)
  par(.pardefault)
  grDevices::palette("default")
}
