# Internal functions for performing dimension reduction and/or standardization, as part of the clusterability R package.

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

#' PCA Methods ---------------------------------------------------

#' Because different machines or implementations of PCA can yield
#' differently signed rotation matrices, and thus scores,
#' we multiply all scores by -1 if the first loading is negative.
#' This should ensure consistent results across machines and between SAS and R implementations,
#' assuming the variables are ordered the same.


#' Compute and return the scores for the first principal component in PCA.
perform_pca <- function(x, center, scale) {
  pcaresult <- stats::prcomp(x, center = center, scale. = scale, retx = TRUE)

  if (pcaresult$rotation[1, 1] < 0) {
    return(-1 * pcaresult$x[, 1])
  } else {
    return(pcaresult$x[, 1])
  }
}

#' Compute and return the scores for the first sparse principal component using the sparsepca implementation.
perform_spca_sparsepca <- function(x, center, scale, alpha, beta) {
  spcaresult <- sparsepca::spca(x, k = 1, alpha = alpha, beta = beta, center = center, scale = scale, verbose = 0)

  if (spcaresult$loadings[1, 1] < 0) {
    return(-1 * spcaresult$scores[, 1])
  } else {
    return(spcaresult$scores[, 1])
  }
}

#' Compute and return the scores for the first sparse principal component using the ELASTICNET implementation.
perform_spca_elasticnet <- function(x, para, lambda) {
  # elasticnet spca uses its own rules for centering and scaling.
  sparse <- "penalty"
  spcaresult <- elasticnet::spca(x, 1, para,
    type = "predictor",
    sparse = sparse, use.corr = FALSE, lambda = lambda,
    max.iter = 200, trace = FALSE, eps.conv = 1e-3
  )

  scores <- x %*% spcaresult$loadings
  if (spcaresult$loadings[1, 1] < 0) {
    return(-1 * scores[, 1])
  } else {
    return(scores[, 1])
  }
}

#' Pairwise distances ---------------------------------------------------

#' Computes pairwise distances and returns a vector of distances.
compute_pairwise_distances <- function(x, method) {
  ### Supported by default in dist() function ###
  # "minkowski(p)" = Minkowski metric, p
  # "euclidean" = Euclidean
  # "maximum" = maximum
  # "manhatten" = manhatten
  # "canberra" = canberra
  # "binary" = binary

  ### Custom metrics ###
  # http://documentation.sas.com/?cdcId=pgmsascdc&cdcVersion=9.4_3.3&docsetId=statug&docsetTarget=statug_distance_details01.htm&locale=en
  # "corr" = correlation
  # "sqeuc" = square euclidean
  # "sqcorr" = squared correlation metric
  # "cov" = covariance metric

  # Check if it's supported by built-in dist()
  distmethods <- c("euclidean", "maximum", "manhattan", "canberra", "binary")
  is_dist_default <- method %in% distmethods

  # Check if it's a Minkowski metric
  minkowski_regex <- "minkowski\\([[:digit:]]*[.]?[[:digit:]]*\\)"
  is_minkowski <- grepl(minkowski_regex, method, ignore.case = TRUE)

  if (is_dist_default) {
    distresult <- as.vector(stats::dist(x = x, method = method))
  } else if (is_minkowski) {
    # Strip out the unnecessary parts to get the "p"
    pattern1 <- "minkowski\\("
    pattern2 <- "\\)"

    out1 <- sub(pattern1, "", method)
    minkowski_p <- as.numeric(sub(pattern2, "", out1))

    distresult <- as.vector(stats::dist(x = x, method = "minkowski", p = minkowski_p))
  } else {
    distresult <- switch(method,
      "sqeuc" = get_distance_squared_euclidean(x),
      "corr" = get_distance_correlation(x),
      "cov" = get_distance_covariance(x),
      "sqcorr" = get_distance_squared_correlation(x)
    )
  }

  return(distresult)
}

#' Returns only complete cases from the original data set.
get_complete_cases <- function(x) {
  ccind <- stats::complete.cases(x)
  return(x[ccind, ])
}

#' Returns the number of rows that are not complete cases.
count_missing_rows <- function(x) {
  total_rows <- NROW(x)
  complete_rows <- NROW(x[stats::complete.cases(x), ])
  return((total_rows - complete_rows))
}

#' Returns a distance matrix using the correlation metric.
get_distance_correlation <- function(x) {
  # Matches SAS
  # Validation handled in validate_metric(). Cannot have 1-dimensional data
  numerator <- ((x - rowMeans(x)) %*% t(x - rowMeans(x)))
  denominator <- rowSums((x - rowMeans(x))^2) %*% t(rowSums((x - rowMeans(x))^2))

  if (any(denominator == 0)) {
    stop("One or more constant observations were found in the data set. When using distance_metric = 'corr', the data must not contain constant observations. ")
  } else {
    final <- numerator / sqrt(denominator)
    # This is a vector
    return(final[lower.tri(final)])
  }
}

#' Returns a distance matrix using the covariance metric.
get_distance_covariance <- function(x) {
  # Matches SAS
  # Validation handled in validate_metric(). Cannot have df = 0
  # This is a similarity metric, not a distance metric.
  df <- dim(x)[2] - 1
  fullmat <- 1 / df * ((x - rowMeans(x)) %*% t(x - rowMeans(x)))
  return(fullmat[lower.tri(fullmat)])
}

#' Returns a distance matrix using the squared euclidean metric.
get_distance_squared_euclidean <- function(x) {
  # Matches SAS
  return(as.vector(stats::dist(x, method = "euclidean"))^2)
}

#' Returns a distance matrix using the squared correlation metric.
get_distance_squared_correlation <- function(x) {
  # Matches SAS
  # Validation handled in validate_metric(). Cannot have 1-dimensional data
  return(get_distance_correlation(x)^2)
}

#' Returns the lower triangular portion of a distance matrix.
get_lower_triangle <- function(x) {
  if (NROW(x) != NCOL(x)) {
    stop("When using the 'is_dist_matrix' argument, the 'data' argument must be a square matrix.")
  }
  return(x[lower.tri(x)])
}
