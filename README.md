Model based clustering for mixed data: clustMD
================
Damien McParland
March 22, 2017

This R package allows the user to perform model based clustering of mixed data (i.e. data that consist of continuous, binary, ordinal or nominal variables) using a parsimonious mixture of latent Gaussian variable models.

This model based clustering approach assumes that underlying the observed categorical response is a latent continuous variable. A finite mixture model is used to identify sub populations or clusters within the larger population.

Installation
------------

The clustMD package can be easily installed in R as follows.

``` r
  install.packages("clustMD")
```

The `Byar` data set that is used in the examples is included in the package. This data set contains information on 475 prostate cancer patients. Measurements taken on these patients consist of continuous, binary, ordinal and nominal variables.

Functions
---------

### `clustMD()`

To use clustMD to cluster the Byar data set you may run the following code. The code consists of some simple pre-processing steps followed by the correct usage of the `clustMD()` function.

``` r
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

    res <- clustMD(X=Y, G=3, CnsIndx=8, OrdIndx=11, Nnorms=20000, 
    MaxIter=500, model="EVI", store.params=FALSE, scale=TRUE, 
    startCL="kmeans")
```

The `clustMD()` function outputs an S3 object of class `clustMD`. Basic S3 methods are included in the package also. The functions available are

-   `print.clustMD()`
-   `summary.clustMD()`
-   `plot.clustMD()`

The `plot.clustMD()` function produces a number of useful summary plots of the `clustMD` object.

### `clustMDparallel()`

Another function is available to run multiple models in parallel called `clustMDparallel()`. This function takes a range of possible values for the number of clusters as a vector. It also takes a character vector as an input that specifies which of the covariance models are to be fitted.

``` r
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

  res <- clustMDparallel(X=Y, G=1:3, CnsIndx=8, OrdIndx=11, Nnorms=20000, 
  MaxIter=500, models=c("EVI", "EII", "VII"), store.params=FALSE, 
  scale=TRUE, startCL="kmeans")
```

The `clustMDparallel()` function outputs an S3 object of class `clustMDparallel`. Some S3 methods are also available for this class:

-   `print.clustMDparallel()`
-   `summary.clustMDparallel()`
-   `plot.clustMDparallel()`

The `plot.clustMDparallel()` function outputs the same plots as `plot.clustMD()` but for the optimal model according to the approximated BIC criterion. An additional plot is also included that illustrated the approximated BIC values for the fitted models.

<!-- README.md is generated from README.Rmd. Please edit that file -->
