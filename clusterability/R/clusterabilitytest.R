#' Perform a test of clusterability
# Copyright (C) 2026  Zachariah Neville, Naomi Brownstein, Andreas Adolfsson, Margareta Ackerman

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#'
#'
#' @description Performs tests for clusterability of a data set and returns results in a clusterability object. Can do data reduction via PCA, sparse PCA, or pairwise distances and standardize data prior to performing the test.
#' @param data the data set to be used in the test. Must contain only numeric data.
#' @param test the test to be performed. Either \code{"dip"} or \code{"silverman"}. See 'Details' section below for how to pick a test.
#' @param reduction any dimension reduction that is to be performed.
#' \itemize{
#' \item{\code{"none"} performs no dimension reduction.}
#' \item{\code{"pca"} uses the scores from the first principal component.}
#' \item{\code{"spca"} uses the scores from the first (sparse) principal component.}
#' \item{\code{"distance"} computes pairwise distances (using \code{distance_metric} as the metric).}
#' }
#' For multivariate \code{data}, dimension reduction is required.
#' @param distance_metric if applicable, the metric to be used in computing pairwise distances.
#'
#' The \code{"euclidean"} (default), \code{"maximum", "manhattan", "canberra", "binary"} choices work the same as in \code{\link{dist}}. The Minkowski metric is available by providing \code{"minkowski(p)"}.
#'
#' Additional choices are:
#' \itemize{
#' \item{\code{"sqeuc"}: squared Euclidean distances.}
#' \item{\code{"cov"}: covariance similarity coefficient,}
#' \item{\code{"corr"}: correlation similarity coefficient}
#' \item{\code{"sqcorr"}: squared correlation similarity coefficient.}
#' }
#'
#' CAUTION: Not all of these have been tested, but instead are provided to potentially be useful. If in doubt, use the default \code{"euclidean"}.
#' @param distance_standardize how the variables should be standardized, if at all.
#' \itemize{
#' \item{\code{"none"}: no standardization is performed}
#' \item{\code{"std"} (default) each variable standardized to have mean 0 and standard deviation 1}
#' \item{\code{"mean"}: each variable standardized to have mean 0 (standard deviation is unchanged)}
#' \item{\code{"median"}: each variable standardized to have median 0 (standard deviation is unchanged)}
#' }
#' @param pca_center if applicable, a logical value indicating whether the variables should be shifted to be zero centered (see \code{\link{prcomp}} for more details). Default is \code{TRUE}.
#' @param pca_scale if applicable, a logical value indicating whether the variables should be scaled to have unit variance before the analysis takes place (see \code{\link{prcomp}} for details). Default is \code{TRUE}.
#' @param spca_method if applicable, the sparse PCA method to use. Either \code{"EN"} for the elasticnet implementation or \code{"VP"} for the variable projection implementation.
#' @param spca_EN_para if applicable, a positive number used as a penalty parameter. Default is \code{0.01}.
#' @param spca_EN_lambda if applicable, the quadratic penalty parameter. Default is \code{1e-6}.
#' @param spca_VP_center if applicable, a logical value indicating whether the variables should be shifted to be zero centered. Default is \code{TRUE}.
#' @param spca_VP_scale if applicable, a logical value indicating whether the variables should be scaled to have unit variance. Default is \code{TRUE}.
#' @param spca_VP_alpha if applicable, the sparsity controlling parameter. Default is \code{1e-3}.
#' @param spca_VP_beta if applicable, the amount of ridge shrinkage to apply in order to improve conditioning. Default is \code{1e-3}.
#' @param is_dist_matrix a logical value indicating whether the \code{data} argument is a distance matrix. If \code{TRUE} then the lower triangular portion of \code{data} will be extracted and be used in the multimodality test.
#' @param completecase a logical value indicating whether a complete case analysis should be performed. For both tests, missing data must be removed before the test can be performed. This can be done manually by the user or by setting \code{completecase = TRUE}.
#' @param d_simulatepvalue for Dip Test, a logical value indicating whether \eqn{p}~values should be obtained via Monte Carlo simulation (see \code{\link[diptest]{dip.test}} for details).
#' @param d_reps for Dip Test, a positive integer. The number of replicates used in Monte Carlo simulation. Only used if \code{d_simulatepvalue} is \code{TRUE}.
#' @param s_m for Silverman Test, a positive integer. The number of bootstrap replicates used in the test. Default is \code{999}.
#' @param s_adjust for Silverman Test, a logical value indicating whether p-values are adjusted using work by Hall and York.
#' @param s_digits for Silverman Test, a positive integer indicating the number of digits to round the p value. Default is \code{6} and is only used when \code{s_adjust = TRUE}.
#' @param s_setseed for Silverman Test, an integer used to set the seed of the random number generator. If the default value of \code{NULL} is used, then no seed will be set.
#' @param s_outseed for Silverman Test, a logical value indicating whether to return the state of the random number generator as part of the output. This is used in limited cases for troubleshooting, so the default is \code{FALSE}.
#' @return \code{clusterabilitytest} returns a \code{clusterability} object containing information on the test performed and results. Can be printed using the \code{\link{print.clusterability}} function.
#' @seealso \code{\link{print.clusterability}}
#'
#'
#' @examples
#' ### Quick start ###
#' # Load data and remove Species
#' data(iris)
#' iris_num <- iris[, -5]
#' plot(iris_num)
#'
#' # Run test using default options
#' clust_result <- clusterabilitytest(iris_num, "dip")
#'
#' # Print results
#' print(clust_result)
#'
#' \donttest{
#' ### Longer Example: Specifying Parameters ###
#' # Load data and plot to visualize
#' data(normals2)
#' plot(normals2)
#'
#' # Using Silverman's test, pairwise distances to reduce dimension,
#' # 1,000 bootstrap replicates, with an RNG seed of 12345
#' clust_result2 <- clusterabilitytest(normals2, "silverman",
#'   reduction = "distance",
#'   s_m = 1000, s_setseed = 12345
#' )
#'
#' # Print result
#' print(clust_result2)
#' }
#'
#' @references Hall, P. and York, M., 2001. On the calibration of Silverman's test for multimodality. Statistica Sinica, pp.515-536.
#'
#' Silverman, B.W., 1981. Using kernel density estimates to investigate multimodality. Journal of the Royal Statistical Society. Series B (Methodological), pp.97-99.
#'
#' Martin Maechler (2016). diptest: Hartigan's Dip Test Statistic for Unimodality - Corrected. R package version 0.75-7. https://CRAN.R-project.org/package=diptest
#'
#' Schwaiger F, Holzmann H. Package which implements the silvermantest; 2013. Available from: https://www.mathematik.uni-marburg.de/stochastik/R packages/.
#'
#' @export
clusterabilitytest <- function(data, test, reduction = "pca", distance_metric = "euclidean",
                               distance_standardize = "std", pca_center = TRUE, pca_scale = TRUE,
                               spca_method = "EN",
                               spca_EN_para = 0.01, spca_EN_lambda = 1e-6,
                               spca_VP_center = TRUE, spca_VP_scale = TRUE,
                               spca_VP_alpha = 1e-3, spca_VP_beta = 1e-3,
                               is_dist_matrix = FALSE, completecase = FALSE, d_simulatepvalue = FALSE,
                               d_reps = 2000, s_m = 999, s_adjust = TRUE, s_digits = 6, s_setseed = NULL, s_outseed = FALSE) {
  # Get name of data set
  dataname <- deparse(substitute(data))

  # Validate data
  validate_data(data, dataname)

  # Gather info about the original data set for printing later
  data <- as.matrix(data)
  nobs <- NROW(data)
  nvar <- NCOL(data)
  nmiss <- count_missing_rows(data)

  # Validate parameters. Necessary before doing any work.
  # Just validate everything every time
  test <- validate_test(test)
  reduction <- validate_reduction(reduction)
  pca_center <- validate_pca_center(pca_center)
  pca_scale <- validate_pca_scale(pca_scale)
  spca_method <- validate_spca_method(spca_method)
  spca_VP_center <- validate_spca_VP_center(spca_VP_center)
  spca_VP_scale <- validate_spca_VP_scale(spca_VP_scale)
  # todo check vp_alpha vp_beta en_para en_lambda

  is_dist_matrix <- validate_isdistmatrix(is_dist_matrix, reduction, data)
  distance_metric <- validate_metric(distance_metric, data)
  distance_standardize <- validate_standardize(distance_standardize)
  completecase <- validate_completecase(completecase)
  d_simulatepvalue <- validate_dsimulatepvalue(d_simulatepvalue)
  d_reps <- validate_dreps(d_reps)
  # s_k <- validate_sk(s_k)
  s_k <- 1
  s_m <- validate_sm(s_m)
  s_adjust <- validate_sadjust(s_adjust)
  s_digits <- validate_sdigits(s_digits)
  s_setseed <- validate_ssetseed(s_setseed)
  s_outseed <- validate_soutseed(s_outseed)

  # Get complete cases
  if (completecase) {
    data <- get_complete_cases(data)
  } else if (nmiss > 0) {
    stop("Missing data was found in the data set and the 'completecase' parameter was not set to TRUE. The Dip and Silverman tests cannot be performed if there is missing data.")
  }

  # Perform standardization if using pairwise distances. Even when "NONE" this works.
  if (identical(reduction, "DISTANCE")) {
    data <- standardize_data(data, distance_standardize)
  }

  # Perform dimension reduction if requested
  if (identical(reduction, "PCA")) {
    data <- perform_pca(data, pca_center, pca_scale)
  } else if (identical(reduction, "SPCA")) {
    if (identical(spca_method, "VP")) {
      data <- perform_spca_sparsepca(
        data, spca_VP_center,
        spca_VP_scale, spca_VP_alpha, spca_VP_beta
      )
    } else if (identical(spca_method, "EN")) {
      data <- perform_spca_elasticnet(data, spca_EN_para, spca_EN_lambda)
    } else {
      stop("Supported SPCA methods are: VP and EN ")
    }
  } else if (identical(reduction, "DISTANCE")) {
    data <- compute_pairwise_distances(data, distance_metric)
  } else if (is_dist_matrix) {
    data <- get_lower_triangle(data)
  }

  # Run the test
  if (identical(test, "DIP")) {
    dipresult <- perform_dip_test(data, d_simulatepvalue, d_reps)
  } else {
    silvresult <- perform_silverman_test(as.matrix(data), s_k, s_m, s_adjust, s_digits, s_setseed)
  }


  # Prepare results for output
  clustobj <- list()
  arglist <- list(data = dataname, test = test, reduction = reduction, completecase = completecase)
  datainfo <- list(numobs = nobs, missingobs = nmiss, numvar = nvar)

  if (identical(reduction, "DISTANCE")) {
    arglist$distance_metric <- distance_metric
    arglist$distance_standardize <- distance_standardize
  } else if (identical(reduction, "PCA")) {
    arglist$pca_center <- pca_center
    arglist$pca_scale <- pca_scale
  } else if (identical(reduction, "SPCA")) {
    arglist$spca_method <- spca_method

    arglist$spca_VP_center <- spca_VP_center
    arglist$spca_VP_scale <- spca_VP_scale
    arglist$spca_VP_alpha <- spca_VP_alpha
    arglist$spca_VP_beta <- spca_VP_beta
    arglist$spca_EN_para <- spca_EN_para
    arglist$spca_EN_lambda <- spca_EN_lambda
  }

  if (identical(test, "DIP")) {
    arglist$simulatepvalues <- d_simulatepvalue

    if (d_simulatepvalue) {
      arglist$reps <- d_reps
    }

    clustobj$nullhypothesis <- "Number of modes = 1"
    clustobj$alternativehypothesis <- "Number of modes > 1"
    clustobj$dipstatistic <- unname(dipresult$statistic)
    clustobj$pvalue <- dipresult$p.value
  } else {
    clustobj$nullhypothesis <- paste("Number of modes <= ", s_k, sep = "")
    clustobj$alternativehypothesis <- paste("Number of modes > ", s_k, sep = "")
    clustobj$critbw <- silvresult$hcrit
    clustobj$k <- silvresult$k
    clustobj$pvalue <- silvresult$p_value

    arglist$adjustpval <- s_adjust
    arglist$reps <- s_m

    if (s_adjust & s_k == 1) {
      arglist$digits <- s_digits
    }

    if (!is.null(s_setseed)) {
      arglist$rngseed <- s_setseed
    }

    if (s_outseed) {
      arglist$rngstate <- silvresult$saved_seed
    }
  }

  clustobj$arglist <- arglist
  clustobj$datainfo <- datainfo

  return(clusterability(clustobj))
}


#' Print a clusterability object
#' @description Print function to display results from a clusterability test.
#' @param x An object of class 'clusterability'
#' @param ... Not used
#'
#'
#' @seealso \code{\link{clusterabilitytest}}
#' @export
print.clusterability <- function(x, ...) {
  # The ... is necessary because the generic print() has it

  correct_class <- tryCatch(identical(class(x), "clusterability"),
    error = function(msg) stop("This function is only meant to print 'clusterability' objects")
  )

  if (!correct_class) {
    stop("This function is only meant to print 'clusterability' objects")
  }

  # Basic info
  cat("----------------------\n")
  cat("Clusterability Test\n")
  cat("----------------------\n\n")

  cat(paste("Data set name: ", x$arglist$data, "\n", sep = ""))
  cat(paste("Your data set has", x$datainfo$numobs, "observation(s) and", x$datainfo$numvar, "variable(s).\n", sep = " "))

  if (x$datainfo$missingobs > 0) {
    cat(paste("WARNING: Of these,", x$datainfo$missingobs, "observation(s) is/are missing at least one variable value.\n\n", sep = " "))
  } else {
    cat("There were no missing values. Your data set is complete.\n\n")
  }


  # Reduction method
  if (identical(x$arglist$reduction, "DISTANCE")) {
    cat("Data Reduced Using: Pairwise Distances\n")
    cat(paste("Distance Metric:", x$arglist$distance_metric, "\n", sep = " "))
    if (!identical(x$arglist$distance_standardize, "NONE")) {
      cat(paste("Standardization Method Used:", x$arglist$distance_standardize, "\n", sep = " "))

      std_description <- switch(toupper(x$arglist$distance_standardize),
        "STD" = "Data is standardized so that each variable has mean 0 and standard deviation 1.",
        "MEAN" = "Each variable is centered around its mean. In other words, the mean - after standardization - is 0, and the standard deviation is unchanged.",
        "MEDIAN" = "Each variable is centered around its median. The median - after standardization - is 0, and the standard deviation is unchanged."
      )

      cat(paste("\t", std_description, sep = " ")) # description
    }
    cat("\n")
  } else if (identical(toupper(x$arglist$reduction), "PCA")) {
    cat("Data Reduced Using: PCA\n")
  } else if (identical(toupper(x$arglist$reduction), "SPCA")) {
    cat("Data Reduced Using: SPCA\n")
  }

  # Test name and results
  if (identical(x$arglist$test, "DIP")) {
    cat("\n-----------------------------------------\n")
    cat("Results: Dip Test of Unimodality\n")
    cat("-----------------------------------------\n\n")
    cat("Null Hypothesis: number of modes = 1\n")
    cat("Alternative Hypothesis: number of modes > 1\n")
    cat(paste("p-value:", x$pvalue, "\n", sep = " "))
    cat(paste("Dip statistic:", x$dipstatistic, "\n\n", sep = " "))
  } else if (identical(x$arglist$test, "SILVERMAN")) {
    cat("\n----------------------------------------------\n")
    cat("Results: Silverman's Critical Bandwidth Test\n")
    cat("----------------------------------------------\n\n")
    cat(paste("Null Hypothesis: number of modes <=", x$k, "\n", sep = " "))
    cat(paste("Alternative Hypothesis: number of modes >", x$k, "\n", sep = " "))
    cat(paste("p-value:", x$pvalue, "\n", sep = " "))
    cat(paste("Critical bandwidth:", x$critbw, "\n\n", sep = " "))
  }

  # Test options used
  cat("---------------------\n")
  cat("Test Options Used\n")
  cat("---------------------\n\n")

  if (identical(x$arglist$test, "DIP")) {
    if (x$arglist$simulatepvalues) {
      cat(paste("p values obtained from a Monte Carlo simulation with", x$arglist$reps, "replicates\n", sep = " "))
      nooptions <- FALSE
    } else {
      cat("Default values for the optional parameters were used. To learn more about customizing the behavior of the clusterabilitytest, please see the R documentation.")
    }
  } else if (identical(x$arglist$test, "SILVERMAN")) {
    if (!is.null(x$arglist$rngseed)) {
      cat(paste("Seed set in R to:", x$arglist$rngseed, "\n", sep = " "))
    }

    cat(paste("p value based on", x$arglist$reps, "bootstrap replicates\n", sep = " "))

    if (x$arglist$adjustpval) {
      cat("Adjusted p-values based on the adjustment in Hall and York (2001)\n")
    }

    if (x$arglist$adjustpval && x$k == 1) {
      cat(paste("p value rounded to:", x$arglist$digits, "digits\n", sep = " "))
    }
  }
}


# Constructor for clusterability objects
# This is only called from within clusterabilitytest() after all parameters were validated
clusterability <- function(paramlist) {
  # Return a new object
  return(structure(paramlist, class = "clusterability"))
}
