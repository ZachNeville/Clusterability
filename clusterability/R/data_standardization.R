# Functions to standardize a data set.
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


# Standardize a data set and return the standardized data
standardize_data <- function(x, method) {
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
