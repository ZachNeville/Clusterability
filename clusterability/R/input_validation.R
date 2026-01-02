#' Internal functions for parameter checking, as part of the clusterability R package
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

#' Validates the data provided
#' @noRd
validate_data <- function(d, data_name) {
  data_is_null <- tryCatch(
    isTRUE(is.null(d)),
    error = function(e) return(TRUE)
  )

  if (data_is_null) {
    stop(paste("The dataset", data_name, "is NULL or is invalid."))
  }

  if (NROW(d) == 0) {
    stop(paste("The dataset", data_name, "has 0 rows."))
  }

  if (NCOL(d) == 0) {
    stop(paste("The dataset", data_name, "has 0 columns."))
  }
}


#' Validates the 'test' parameter and returns an uppercase version of it
#' @noRd
validate_test <- function(test) {
  valid_tests <- c("DIP", "SILVERMAN")

  test_is_valid <- tryCatch(
    isTRUE((!is.null(test) && (toupper(test) %in% valid_tests))),
    error = function(e) return(FALSE)
  )

  if (!test_is_valid) {
    stop("Invalid `test` argument was specified. The `test` parameter must be either \"dip\" or \"silverman\"")
  } else {
    toupper(test)
  }
}

#' Validates the 'distance_metric' parameter and returns a lowercase version of it
#' @noRd
validate_metric <- function(metric, x) {
  metric_is_null <- tryCatch(
    isTRUE(is.null(metric)),
    error = function(e) return(TRUE)
  )

  if (metric_is_null) {
    warning("Invalid distance metric was entered. The default metric (\"euclidean\") will be used. Please see the help file for a list of valid metrics.")
    return("euclidean")
  } else {
    metric_lowercase <- tolower(metric)

    valid_metrics <- c("euclidean", "maximum", "manhattan", "canberra", "binary", "sqeuc", "sqcorr", "corr", "cov")

    # Check for a valid minkowski metric. Form must be "minkowski(p)", where p is a positive numeric value (not necessarily integer)
    minkowski_regex <- "minkowski\\([[:digit:]]*[.]?[[:digit:]]*\\)"
    is_minkowski_metric <- grepl(minkowski_regex, metric_lowercase, ignore.case = TRUE)

    # If it's one of the valid metrics, return the lowercase version
    if (isTRUE(metric_lowercase %in% valid_metrics)) {
      if (isTRUE(metric_lowercase %in% c("cov", "corr", "sqcorr") && identical(as.double(NCOL(x)), 1))) {
        stop('The "cov", "corr", and "sqcorr" metrics are not available for 1-dimensional data.')
      } else {
        return(metric_lowercase)
      }
    } else if (is_minkowski_metric) {
      # If it's a minkowski metric, then check the "p" value
      pattern1 <- "minkowski\\("
      pattern2 <- "\\)"

      out1 <- sub(pattern1, "", metric_lowercase)
      minkowski_p <- as.numeric(sub(pattern2, "", out1))

      # "p" must be positive numeric
      if (minkowski_p > 0) {
        return(metric_lowercase)
      } else {
        warning("Invalid value for p was entered when using the Minkowski metric. p must be a positive number. The default value of 2 will be used.")
        return("minkowski(2)")
      }
    } else {
      warning("Invalid distance metric was entered. The default metric (\"euclidean\") will be used. Please see the help file for a list of valid metrics.")
      return("euclidean")
    }
  }
}

#' Validates the 'reduction' argument and returns an uppercase version of it
#' @noRd
validate_reduction <- function(reduction) {
  valid_reductions <- c("PCA", "SPCA", "DISTANCE", "NONE")
  reduction_is_valid <- tryCatch(
    isTRUE((!is.null(reduction) && (toupper(reduction) %in% valid_reductions))),
    error = function(e) return(FALSE)
  )

  if (!reduction_is_valid) {
    warning("Invalid reduction method was used. No reduction was performed. The `reduction` argument must be \"PCA\", \"SPCA\", \"DISTANCE\", or \"NONE\"")
    return("NONE")
  } else {
    toupper(reduction)
  }
}

#' Validates the 'is_dist_matrix' argument and returns its value
#' @noRd
validate_isdistmatrix <- function(is_dist_matrix, reduction, data) {
  if (validate_boolean(is_dist_matrix, "is_dist_matrix", FALSE)) {
    if (!identical(reduction, "NONE")) {
      stop("When providing a value of TRUE for the `is_dist_matrix`` argument, the `reduction` argument must be \"NONE\".")
    } else if (!isSymmetric(as.matrix(data), check.attributes = FALSE)) {
      stop("When providing a value of TRUE for the `is_dist_matrix` argument, the dataset must be a symmetric matrix.")
    } else {
      return(TRUE)
    }
  } else {
    FALSE
  }
}


#' Validates the 'distance_standardize' argument and returns an uppercase version of it
#' @noRd
validate_standardize <- function(standard) {
  valid_standardize_values <- c("STD", "NONE", "MEAN", "MEDIAN")
  is_valid <- tryCatch(
    isTRUE((!is.null(standard) && (toupper(standard) %in% valid_standardize_values))),
    error = function(e) return(FALSE)
  )

  if (!is_valid) {
    warning("Invalid standardization technique was used. No standardization was performed. Please see documentation for valid techniques.")
    return("NONE")
  } else {
    toupper(standard)
  }
}

#' Validate that pca_center is a logical type
#' @noRd
validate_pca_center <- function(center) {
  validate_boolean(center, "pca_center", TRUE)
}

#' Validate that pca_scale is a logical type
#' @noRd
validate_pca_scale <- function(scale) {
  validate_boolean(scale, "pca_scale", TRUE)
}


#' Validates the 'spca_method' argument and returns an uppercase version of it
#' @noRd
validate_spca_method <- function(spca_method) {
  valid_spca_methods <- c("EN", "VP")
  is_valid <- tryCatch(
    isTRUE((!is.null(spca_method) && (toupper(spca_method) %in% valid_spca_methods))),
    error = function(e) return(FALSE)
  )

  if (!is_valid) {
    stop("Invalid spca method was used. The `spca_method` argument must be \"EN\" or \"VP\"")
  } else {
    toupper(spca_method)
  }
}

#' Validate that spca_VP_center is a logical type
#' @noRd
validate_spca_VP_center <- function(center) {
  validate_boolean(center, "spca_VP_center", TRUE)
}

#' Validate that spca_VP_scale is a logical type
#' @noRd
validate_spca_VP_scale <- function(scale) {
  validate_boolean(scale, "spca_VP_scale", TRUE)
}

#' Validates the spca_VP_alpha parameter
#' @noRd
validate_spca_VP_alpha <- function(spca_VP_alpha){
  is_valid <- tryCatch(
    isTRUE((!is.null(spca_VP_alpha) && is.numeric(spca_VP_alpha))),
    error = function(e) return(FALSE)
  )

  if (!is_valid) {
    warning("The value of `spca_VP_alpha` must be numeric. The default value of 1e-3 will be used.")
    return(1e-3)
  } else {
    spca_VP_alpha
  }
}

#' Validates the spca_VP_beta parameter
#' @noRd
validate_spca_VP_beta <- function(spca_VP_beta){
  is_valid <- tryCatch(
    isTRUE((!is.null(spca_VP_beta) && is.numeric(spca_VP_beta))),
    error = function(e) return(FALSE)
  )

  if (!is_valid) {
    warning("The value of `spca_VP_beta` must be numeric. The default value of 1e-3 will be used.")
    return(1e-3)
  } else {
    spca_VP_beta
  }
}

#' Validates the spca_EN_para parameter
#' @noRd
validate_spca_EN_para <- function(spca_EN_para){
  is_valid <- tryCatch(
    isTRUE((!is.null(spca_EN_para) && is.numeric(spca_EN_para) && spca_EN_para > 0)),
    error = function(e) return(FALSE)
  )

  if (!is_valid) {
    warning("The value of `spca_EN_para` must be a positive numeric value. The default value of 0.01 will be used.")
    return(0.01)
  } else {
    spca_EN_para
  }
}

#' Validates the spca_EN_lambda parameter
#' @noRd
validate_spca_EN_lambda <- function(spca_EN_lambda){
  is_valid <- tryCatch(
    isTRUE((!is.null(spca_EN_lambda) && is.numeric(spca_EN_lambda))),
    error = function(e) return(FALSE)
  )

  if (!is_valid) {
    warning("The value of `spca_EN_lambda` must be numeric. The default value of 1e-6 will be used.")
    return(1e-6)
  } else {
    spca_EN_lambda
  }
}

#' Validate that completecase is a logical type
#' @noRd
validate_completecase <- function(completecase) {
  validate_boolean(completecase, "completecase", FALSE)
}

#' Validate that d_simulatepvalue is a logical type
#' @noRd
validate_dsimulatepvalue <- function(d_simulatepvalue) {
  validate_boolean(d_simulatepvalue, "d_simulatepvalue", FALSE)
}

#' Validate that s_adjust is a logical type
#' @noRd
validate_sadjust <- function(s_adjust) {
  validate_boolean(s_adjust, "s_adjust", TRUE)
}

#' Validate that s_outseed is a logical type
#' @noRd
validate_soutseed <- function(s_outseed) {
  validate_boolean(s_outseed, "s_outseed", FALSE)
}

#' Validate d_reps parameter
#' @noRd
validate_dreps <- function(d_reps) {
  is_valid <- tryCatch(
    isTRUE((!is.null(d_reps) && is.numeric(d_reps) && (d_reps >= 1) && ((d_reps %% 1) == 0))),
    error = function(e) return(FALSE)
  )

  if (!is_valid) {
    warning("The value of `d_reps` must be a positive integer. The default value of 2000 will be used.")
    return(2000)
  } else {
    d_reps
  }
}

#' Validate s_m parameter
#' @noRd
validate_sm <- function(s_m) {
  is_valid <- tryCatch(
    isTRUE((!is.null(s_m) && is.numeric(s_m) && (s_m >= 1) && ((s_m %% 1) == 0))),
    error = function(e) return(FALSE)
  )

  if (!is_valid) {
    warning("The value of `s_m` must be a positive integer. The default value of 999 will be used.")
    return(999)
  } else {
    s_m
  }
}

#' Validate s_digits parameter
#' @noRd
validate_sdigits <- function(s_digits) {
  is_valid <- tryCatch(
    isTRUE((!is.null(s_digits) && is.numeric(s_digits) && (s_digits >= 1) && ((s_digits %% 1) == 0))),
    error = function(e) return(FALSE)
  )

  if (!is_valid) {
    warning("The value of `s_digits` must be a positive integer. The default value of 6 will be used.")
    return(6)
  } else {
    s_digits
  }
}

#' Validate the s_setseed parameter
#' @noRd
validate_ssetseed <- function(s_setseed) {
  # Note that NULL is the default so it should still be allowed
  is_integer <- tryCatch(
    isTRUE((is.numeric(s_setseed) && (s_setseed %% 1 == 0))),
    error = function(e) return(FALSE)
  )

  is_null <- tryCatch(
    is.null(s_setseed),
    error = function(e) return(FALSE)
  )

  if (!is_integer && !is_null) {
    warning("The value of `s_setseed` must be an integer. The seed was not set.")
    return(NULL)
  } else {
    s_setseed
  }
}


#' Helper function used to validate a logical (boolean) type
#' @noRd
validate_boolean <- function(var, name, default) {
  is_logical <- tryCatch(
    identical(typeof(var), "logical"),
    error = function(errmsg) return(FALSE)
  )

  if (!is_logical) {
    warning(paste("The `", name, "` argument must be boolean (either TRUE or FALSE). The default value (", default, ") will be used.", sep = ""))
    return(default)
  } else {
    var
  }
}
