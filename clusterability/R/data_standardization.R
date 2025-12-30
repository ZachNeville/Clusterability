# Standardize a data set and return the standardized data
standardizedata <- function(x, method) {
  # NONE = don't do anything
  # STD = mean 0, stdev 1
  # MEAN = subtract mean
  # MEDIAN = subtract median
  # These match with SAS results
  # http://documentation.sas.com/?cdcId=pgmsascdc&cdcVersion=9.4_3.3&docsetId=statug&docsetTarget=statug_stdize_details01.htm&locale=en

  result <- switch(method,
    "STD" = scale(x),
    "NONE" = x,
    "MEAN" = apply(x, 2, function(x) (x - mean(x))),
    "MEDIAN" = apply(x, 2, function(x) (x - stats::median(x)))
  )

  if (any(is.nan(result))) {
    warning("NaN values occurred during standardization. One possible cause is that the data contains a variable which is constant. No standardization was performed.")
    return(x)
  } else {
    return(result)
  }
}
