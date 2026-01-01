# Timings for the clusterability R package

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

##### Initial setup #####
library(clusterability)
library(bench)
ntimes <- 10

test_and_print <- function(data, test, reduction, seed, spca_method = "EN") {
  if (test == "dip") {
    benchmarkResult <- bench::mark(clusterabilitytest(data, "dip", reduction, spca_method = spca_method), min_iterations = ntimes)
  } else if (test == "silverman") {
    benchmarkResult <- bench::mark(clusterabilitytest(data, "silverman", reduction, spca_method = spca_method, distance_standardize = "NONE", s_setseed = seed), min_iterations = ntimes)
  } else {
    stop("Invalid test")
  }

  print(paste("Median Time: ", bench::as_bench_time(benchmarkResult["median"]), ". Iterations: ", benchmarkResult["n_itr"], sep = " "))
}

##### normals1 #####
data(normals1)
normals1 <- normals1[, -3]
test_and_print(normals1, "dip", "pca", NULL)
test_and_print(normals1, "dip", "distance", NULL)
test_and_print(normals1, "silverman", "pca", 123)
test_and_print(normals1, "silverman", "distance", 123)

test_and_print(normals1, "dip", "spca", NULL, "EN")
test_and_print(normals1, "dip", "spca", NULL, "VP")
test_and_print(normals1, "silverman", "spca", 123, "EN")
test_and_print(normals1, "silverman", "spca", 123, "VP")

##### normals2 #####
data(normals2)
normals2 <- normals2[, -3]
test_and_print(normals2, "dip", "pca", NULL)
test_and_print(normals2, "dip", "distance", NULL)
test_and_print(normals2, "silverman", "pca", 123)
test_and_print(normals2, "silverman", "distance", 123)

test_and_print(normals2, "dip", "spca", NULL, "EN")
test_and_print(normals2, "dip", "spca", NULL, "VP")
test_and_print(normals2, "silverman", "spca", 123, "EN")
test_and_print(normals2, "silverman", "spca", 123, "VP")

##### normals3 #####
data(normals3)
normals3 <- normals3[, -3]
test_and_print(normals3, "dip", "pca", NULL)
test_and_print(normals3, "dip", "distance", NULL)
test_and_print(normals3, "silverman", "pca", 123)
test_and_print(normals3, "silverman", "distance", 123)

test_and_print(normals3, "dip", "spca", NULL, "EN")
test_and_print(normals3, "dip", "spca", NULL, "VP")
test_and_print(normals3, "silverman", "spca", 123, "EN")
test_and_print(normals3, "silverman", "spca", 123, "VP")


##### normals4 #####
data(normals4)
normals4 <- normals4[, -4]
test_and_print(normals4, "dip", "pca", NULL)
test_and_print(normals4, "dip", "distance", NULL)
test_and_print(normals4, "silverman", "pca", 123)
test_and_print(normals4, "silverman", "distance", 123)

test_and_print(normals4, "dip", "spca", NULL, "EN")
test_and_print(normals4, "dip", "spca", NULL, "VP")
test_and_print(normals4, "silverman", "spca", 123, "EN")
test_and_print(normals4, "silverman", "spca", 123, "VP")

##### normals5 #####
data(normals5)
normals5 <- normals5[, -4]
test_and_print(normals5, "dip", "pca", NULL)
test_and_print(normals5, "dip", "distance", NULL)
test_and_print(normals5, "silverman", "pca", 123)
test_and_print(normals5, "silverman", "distance", 123)

test_and_print(normals5, "dip", "spca", NULL, "EN")
test_and_print(normals5, "dip", "spca", NULL, "VP")
test_and_print(normals5, "silverman", "spca", 123, "EN")
test_and_print(normals5, "silverman", "spca", 123, "VP")

##### cars #####
data(cars)
test_and_print(cars, "dip", "pca", NULL)
test_and_print(cars, "dip", "distance", NULL)
test_and_print(cars, "silverman", "pca", 123)
test_and_print(cars, "silverman", "distance", 123)

test_and_print(cars, "dip", "spca", NULL, "EN")
test_and_print(cars, "dip", "spca", NULL, "VP")
test_and_print(cars, "silverman", "spca", 123, "EN")
test_and_print(cars, "silverman", "spca", 123, "VP")


##### iris #####
data(iris)
iris_numeric <- iris[, c(1:4)]
test_and_print(iris_numeric, "dip", "pca", NULL)
test_and_print(iris_numeric, "dip", "distance", NULL)
test_and_print(iris_numeric, "silverman", "pca", 123)
test_and_print(iris_numeric, "silverman", "distance", 123)

test_and_print(iris_numeric, "dip", "spca", NULL, "EN")
test_and_print(iris_numeric, "dip", "spca", NULL, "VP")
test_and_print(iris_numeric, "silverman", "spca", 123, "EN")
test_and_print(iris_numeric, "silverman", "spca", 123, "VP")
