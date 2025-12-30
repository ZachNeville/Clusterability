# Code to produce plots for the clusterability package

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


###### Initial Setup ######
library(clusterability)
data(normals1)
data(normals2)
data(normals3)

col1 <- rgb(63, 116, 255, maxColorValue = 255)
col2 <- rgb(225, 4, 43, maxColorValue = 255)
col3 <- rgb(63, 145, 70, maxColorValue = 255)

sym1 <- 21
sym2 <- 23
sym3 <- 24

n1 <- ifelse(normals1$cluster == 1, col1, ifelse(normals1$cluster == 2, col2, col3))
n2 <- ifelse(normals2$cluster == 1, col1, ifelse(normals2$cluster == 2, col2, col3))
n3 <- ifelse(normals3$cluster == 1, col1, ifelse(normals3$cluster == 2, col2, col3))

s1 <- ifelse(normals1$cluster == 1, sym1, ifelse(normals1$cluster == 2, sym2, sym3))
s2 <- ifelse(normals2$cluster == 1, sym1, ifelse(normals2$cluster == 2, sym2, sym3))
s3 <- ifelse(normals3$cluster == 1, sym1, ifelse(normals3$cluster == 2, sym2, sym3))


##### Normals 1 #####
dev.new(width = 300, height = 300, unit = "px") # new window
par(mar = c(5, 5, 2, 2)) # margins
par(xpd = T, mar = par()$mar + c(4.25, 0, 0.75, 0))
plot(normals1[, 1:2], col = "black", pch = s1, cex = 3, cex.main = 3, cex.axis = 1.75, cex.lab = 2.75, main = "normals1", pty = "s", bg = n1)
legend(x = -0.65, y = 0, legend = c("Cluster 1"), col = "black", pch = c(sym1), pt.bg = c(col1), cex = 2, pt.cex = 3, horiz = TRUE)


##### Normals 2 #####
dev.new(width = 300, height = 300, unit = "px")
par(mar = c(5, 5, 2, 2))
par(xpd = T, mar = par()$mar + c(4.25, 0, 0.75, 0))
plot(normals2[, 1:2], col = "black", pch = s2, cex = 3, cex.main = 3, cex.axis = 1.75, cex.lab = 2.75, main = "normals2", pty = "s", bg = n2)
legend(x = -4, y = -6.80, legend = c("Cluster 1", "Cluster 2"), col = "black", pch = c(sym1, sym2), pt.bg = c(col1, col2), cex = 2, pt.cex = 3, horiz = TRUE)


##### Normals 3 #####
dev.new(width = 300, height = 300, unit = "px")
par(mar = c(5, 5, 2, 2))
par(xpd = T, mar = par()$mar + c(4.25, 0, 0.75, 0))
plot(normals3[, 1:2], col = "black", pch = s3, cex = 3, cex.main = 3, cex.axis = 1.75, cex.lab = 2.75, main = "normals3", pty = "s", bg = n3)
legend(x = -2.5, y = -5.75, legend = c("Cluster 1", "Cluster 2", "Cluster 3"), col = "black", pch = c(sym1, sym2, sym3), pt.bg = c(col1, col2, col3), cex = 2, pt.cex = 3, horiz = TRUE)


##### Setup for Normals 4 and 5 #####
library(plotly)
data(normals4)
normals4$markers[which(normals4$cluster == 2)] <- "square"
normals4$markers[which(normals4$cluster == 1)] <- "circle"

data(normals5)
normals5$markers[which(normals5$cluster == 1)] <- "circle"
normals5$markers[which(normals5$cluster == 2)] <- "square"
normals5$markers[which(normals5$cluster == 3)] <- "diamond"

colors2 <- c(col1, col2)
colors3 <- c(col1, col2, col3)

scene <- list(
  xaxis = list(titlefont = list(size = 20), tickfont = list(size = 14)),
  yaxis = list(titlefont = list(size = 20), tickfont = list(size = 14)),
  zaxis = list(titlefont = list(size = 20), tickfont = list(size = 14)),
  camera = list(eye = list(x = 0.8, y = -2.2, z = 0.15))
)


normals4$clustername[which(normals4$cluster == 2)] <- "Cluster 2"
normals4$clustername[which(normals4$cluster == 1)] <- "Cluster 1"


normals5$clustername[which(normals5$cluster == 1)] <- "Cluster 1"
normals5$clustername[which(normals5$cluster == 2)] <- "Cluster 2"
normals5$clustername[which(normals5$cluster == 3)] <- "Cluster 3"


##### Normals 4 #####
plot_ly(data = normals4) %>%
  add_markers(
    x = ~x, y = ~y, z = ~z, type = "scatter3d", mode = "markers", stroke = I("black"), span = I(1.5),
    color = ~clustername, colors = colors2, symbol = ~ I(markers)
  ) %>%
  layout(
    title = list(text = "<b>normals4</b>", yanchor = "bottom", y = 0.80, font = list(size = 24)),
    scene = scene,
    legend = list(font = list(size = 20), orientation = "h", x = 0.2, y = 0.05, borderwidth = 1)
  )


##### Normals 5 #####
plot_ly(data = normals5) %>%
  add_markers(
    x = ~x, y = ~y, z = ~z, type = "scatter3d", mode = "markers", stroke = I("black"), span = I(1.5),
    color = ~clustername, colors = colors3, symbol = ~ I(markers)
  ) %>%
  layout(
    title = list(text = "<b>normals5</b>", yanchor = "bottom", y = 0.85, font = list(size = 24)), scene = scene,
    legend = list(font = list(size = 20), orientation = "h", x = 0.15, y = 0, borderwidth = 1)
  )


##### iris #####
data(iris)

iris <- iris[, -5]
par(mar = c(5, 5, 2, 2))

irispca <- clusterability:::performpca(iris, TRUE, TRUE)

hist(irispca, col = col2, main = "iris", xlab = "Score of First Principal Component", ylab = "Frequency", cex.main = 1.5, cex.axis = 1, cex.lab = 1.25)

irisdist <- as.vector(dist(iris))

hist(irisdist, col = col3, main = "iris", xlab = "Pairwise Distance", ylab = "Frequency", cex.main = 1.5, cex.axis = 1, cex.lab = 1.25)

irisspca <- clusterability:::performspca.elasticnet(as.matrix(iris), 0.01, 1e-6)

hist(irisspca, col = col1, main = "iris", xlab = "Score of First Sparse Principal Component", ylab = "Frequency", cex.main = 1.5, cex.axis = 1, cex.lab = 1.25)


##### cars #####
data(cars)

par(mar = c(5, 5, 2, 2))
plot(cars, col = "black", pch = sym1, cex = 1.5, cex.main = 1.5, cex.axis = 1, cex.lab = 1.25, main = "cars", pty = "s", bg = col1, xlab = "Speed", ylab = "Distance")

carspca <- clusterability:::performpca(cars, TRUE, TRUE)

hist(carspca, col = col1, main = "cars", xlab = "Score of First Principal Component", ylab = "Frequency", cex.main = 1.5, cex.axis = 1, cex.lab = 1.25)

carsdist <- as.vector(dist(cars))

hist(carsdist, col = col2, main = "cars", xlab = "Pairwise Distance", ylab = "Frequency", cex.main = 1.5, cex.axis = 1, cex.lab = 1.25)

carsspca <- clusterability:::performspca.elasticnet(as.matrix(cars), 0.01, 1e-6)

hist(carsspca, col = col3, main = "cars", xlab = "Score of First Sparse Principal Component", ylab = "Frequency", cex.main = 1.5, cex.axis = 1, cex.lab = 1.25)
