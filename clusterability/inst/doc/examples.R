#################################################################################
# Examples file to replicate examples and plots from clusterability package

# Copyright (C) 2025  Zachariah Neville, Naomi Brownstein, Andreas Adolfsson, Margareta Ackerman

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

library(clusterability)

# Normals1 - results
data(normals1)
normals1 <- normals1[, -3]
norm1_dippca <- clusterabilitytest(normals1, "dip")
norm1_dipdist <- clusterabilitytest(normals1, "dip", distance_standardize = "NONE", reduction = "distance")
norm1_silvpca <- clusterabilitytest(normals1, "silverman", s_setseed = 123)
norm1_silvdist <- clusterabilitytest(normals1, "silverman", distance_standardize = "NONE", reduction = "distance", s_setseed = 123)

norm1_dipspcaEN <- clusterabilitytest(normals1, "dip", reduction = "spca", spca_method = "EN")
norm1_dipspcaVP <- clusterabilitytest(normals1, "dip", reduction = "spca", spca_method = "VP")
norm1_silvspcaEN <- clusterabilitytest(normals1, "silverman", s_setseed = 123, reduction = "spca", spca_method = "EN")
norm1_silvspcaVP <- clusterabilitytest(normals1, "silverman", s_setseed = 123, reduction = "spca", spca_method = "VP")


print(norm1_dippca)
print(norm1_dipdist)
print(norm1_silvpca)
print(norm1_silvdist)

print(norm1_dipspcaEN)
print(norm1_dipspcaVP)
print(norm1_silvspcaEN)
print(norm1_silvspcaVP)


# Normals2 - results
data(normals2)
normals2 <- normals2[, -3]
norm2_dippca <- clusterabilitytest(normals2, "dip")
norm2_dipdist <- clusterabilitytest(normals2, "dip", reduction = "distance", distance_standardize = "NONE")
norm2_silvpca <- clusterabilitytest(normals2, "silverman", s_setseed = 123)
norm2_silvdist <- clusterabilitytest(normals2, "silverman", reduction = "distance", distance_standardize = "NONE", s_setseed = 123)

norm2_dipspcaEN <- clusterabilitytest(normals2, "dip", reduction = "spca", spca_method = "EN")
norm2_dipspcaVP <- clusterabilitytest(normals2, "dip", reduction = "spca", spca_method = "VP")
norm2_silvspcaEN <- clusterabilitytest(normals2, "silverman", s_setseed = 123, reduction = "spca", spca_method = "EN")
norm2_silvspcaVP <- clusterabilitytest(normals2, "silverman", s_setseed = 123, reduction = "spca", spca_method = "VP")


print(norm2_dippca)
print(norm2_dipdist)
print(norm2_silvpca)
print(norm2_silvdist)

print(norm2_dipspcaEN)
print(norm2_dipspcaVP)
print(norm2_silvspcaEN)
print(norm2_silvspcaVP)


# Normals3 - results
data(normals3)
normals3 <- normals3[, -3]
norm3_dippca <- clusterabilitytest(normals3, "dip")
norm3_dipdist <- clusterabilitytest(normals3, "dip", reduction = "distance", distance_standardize = "NONE")
norm3_silvpca <- clusterabilitytest(normals3, "silverman", s_setseed = 123)
norm3_silvdist <- clusterabilitytest(normals3, "silverman", reduction = "distance", distance_standardize = "NONE", s_setseed = 123)

norm3_dipspcaEN <- clusterabilitytest(normals3, "dip", reduction = "spca", spca_method = "EN")
norm3_dipspcaVP <- clusterabilitytest(normals3, "dip", reduction = "spca", spca_method = "VP")
norm3_silvspcaEN <- clusterabilitytest(normals3, "silverman", s_setseed = 123, reduction = "spca", spca_method = "EN")
norm3_silvspcaVP <- clusterabilitytest(normals3, "silverman", s_setseed = 123, reduction = "spca", spca_method = "VP")


print(norm3_dippca)
print(norm3_dipdist)
print(norm3_silvpca)
print(norm3_silvdist)

print(norm3_dipspcaEN)
print(norm3_dipspcaVP)
print(norm3_silvspcaEN)
print(norm3_silvspcaVP)


# Normals4 - results
data(normals4)
normals4 <- normals4[, -4]
norm4_dippca <- clusterabilitytest(normals4, "dip")
norm4_dipdist <- clusterabilitytest(normals4, "dip", reduction = "distance", distance_standardize = "NONE")
norm4_silvpca <- clusterabilitytest(normals4, "silverman", s_setseed = 123)
norm4_silvdist <- clusterabilitytest(normals4, "silverman", reduction = "distance", distance_standardize = "NONE", s_setseed = 123)

norm4_dipspcaEN <- clusterabilitytest(normals4, "dip", reduction = "spca", spca_method = "EN")
norm4_dipspcaVP <- clusterabilitytest(normals4, "dip", reduction = "spca", spca_method = "VP")
norm4_silvspcaEN <- clusterabilitytest(normals4, "silverman", s_setseed = 123, reduction = "spca", spca_method = "EN")
norm4_silvspcaVP <- clusterabilitytest(normals4, "silverman", s_setseed = 123, reduction = "spca", spca_method = "VP")


print(norm4_dippca)
print(norm4_dipdist)
print(norm4_silvpca)
print(norm4_silvdist)

print(norm4_dipspcaEN)
print(norm4_dipspcaVP)
print(norm4_silvspcaEN)
print(norm4_silvspcaVP)


# Normals5 - results
data(normals5)
normals5 <- normals5[, -4]
norm5_dippca <- clusterabilitytest(normals5, "dip")
norm5_dipdist <- clusterabilitytest(normals5, "dip", reduction = "distance", distance_standardize = "NONE")
norm5_silvpca <- clusterabilitytest(normals5, "silverman", s_setseed = 123, s_digits = 20)
norm5_silvdist <- clusterabilitytest(normals5, "silverman", reduction = "distance", distance_standardize = "NONE", s_setseed = 123)


norm5_dipspcaEN <- clusterabilitytest(normals5, "dip", reduction = "spca", spca_method = "EN")
norm5_dipspcaVP <- clusterabilitytest(normals5, "dip", reduction = "spca", spca_method = "VP")
norm5_silvspcaEN <- clusterabilitytest(normals5, "silverman", s_setseed = 123, reduction = "spca", spca_method = "EN")
norm5_silvspcaVP <- clusterabilitytest(normals5, "silverman", s_setseed = 123, reduction = "spca", spca_method = "VP")


print(norm5_dippca)
print(norm5_dipdist)
print(norm5_silvpca)
print(norm5_silvdist)

print(norm5_dipspcaEN)
print(norm5_dipspcaVP)
print(norm5_silvspcaEN)
print(norm5_silvspcaVP)


# cars - results
data(cars)

cars_dippca <- clusterabilitytest(cars, "dip")
cars_dipdist <- clusterabilitytest(cars, "dip", reduction = "distance", distance_standardize = "NONE")
cars_silvpca <- clusterabilitytest(cars, "silverman", s_setseed = 123)
cars_silvdist <- clusterabilitytest(cars, "silverman", reduction = "distance", distance_standardize = "NONE", s_setseed = 123)

cars_dipspcaEN <- clusterabilitytest(cars, "dip", reduction = "spca", spca_method = "EN")
cars_dipspcaVP <- clusterabilitytest(cars, "dip", reduction = "spca", spca_method = "VP")
cars_silvspcaEN <- clusterabilitytest(cars, "silverman", s_setseed = 123, reduction = "spca", spca_method = "EN")
cars_silvspcaVP <- clusterabilitytest(cars, "silverman", s_setseed = 123, reduction = "spca", spca_method = "VP")


print(cars_dippca)
print(cars_dipdist)
print(cars_silvpca)
print(cars_silvdist)

print(cars_dipspcaEN)
print(cars_dipspcaVP)
print(cars_silvspcaEN)
print(cars_silvspcaVP)


# iris - results
data(iris)
newiris <- iris[, c(1:4)]
iris_dippca <- clusterabilitytest(newiris, "dip")
iris_dipdist <- clusterabilitytest(newiris, "dip", reduction = "distance", distance_standardize = "NONE")
iris_silvpca <- clusterabilitytest(newiris, "silverman", s_setseed = 123)
iris_silvdist <- clusterabilitytest(newiris, "silverman", reduction = "distance", distance_standardize = "NONE", s_setseed = 123)

iris_dipspcaEN <- clusterabilitytest(newiris, "dip", reduction = "spca", spca_method = "EN")
iris_dipspcaVP <- clusterabilitytest(newiris, "dip", reduction = "spca", spca_method = "VP")
iris_silvspcaEN <- clusterabilitytest(newiris, "silverman", s_setseed = 123, reduction = "spca", spca_method = "EN")
iris_silvspcaVP <- clusterabilitytest(newiris, "silverman", s_setseed = 123, reduction = "spca", spca_method = "VP")

print(iris_dippca)
print(iris_dipdist)
print(iris_silvpca)
print(iris_silvdist)

print(iris_dipspcaEN)
print(iris_dipspcaVP)
print(iris_silvspcaEN)
print(iris_silvspcaVP)
