# Construction des figures
# ************************

# R1+R2 gardées dans myspca1_2
# R1 gardées dans myspca1
# R2 dans myspca2
# males 1 dans myspcam1
# males 2 dans myspcam2

library(maptools)
library(raster)
library(adegenet)

Azrou_shp <- readShapePoly("Azrou_2015_2.shp")

# Plot pour lignée 1 et 2 ensemble
# ********************************

source("sw.value.R")
sc <- as.vector(myspca1_2$li)[,1]
for(i in 1:length(sc)){ 
   if(sc[i] < 0) sc[i] = 0
   if(sc[i] > 0) sc[i] = 1
}
color <- sc
for(i in 1: length(color)){
   if(color[i] == 0) color[i] = "black"
   if(color[i] == 1) color[i] = "white"
}

par(mar = c(0, 0, 0, 0))
plot(NULL, xlim = c(-5.1538, -5.1416), ylim = c(33.426, 33.43385))
plot(Azrou_shp, col = "lightgrey", add = T)
for(i in 1:length(sc))   points(x = spa[i, 2], y = spa[i, 1], pch = 21, bg = color[i])
scalebar(0.1, xy = c(-5.143, 33.4265), type = "bar", divs = 2, below = "kilometres")

# Plot pour la sPCA (par lignée)
# ******************************

par(mar = c(0, 0, 0, 0))
plot(NULL, xlim = c(-5.1538, -5.1416), ylim = c(33.426, 33.43385))
plot(Azrou_shp, col = "lightgrey", add = T)
# Adapter premier argument en fonction de la lignée et du sexe qu'on veut visualiser
# (v. noms des objets en haut de ce script)
colorplot(myspca2, axes = 1:3, add.plot = TRUE, cex = 1.5)

scalebar(0.1, xy = c(-5.143, 33.4265), type = "bar", divs = 2, below = "kilometres")

# Plot pour le clustering
# ***********************

# R1
spa <- read.table("spa_R1_na.txt", h = T, sep = "\t", row.names = 1)
spa <- cbind(spa[, 2], spa[, 1])
sc <- myspca1$li   # Only the selected sPCA axes
spe.dh <- dist(sc)
      # Ward's Minimum Variance Clustering
spe.h.ward <- hclust(spe.dh, method = "ward.D")
color <- c("darkorchid", "firebrick", "black", "chartreuse4", "darksalmon")
shape <- c(21:25)

# R2
spa <- read.table("spa_R2_na.txt", h = T, sep = "\t", row.names = 1)
spa <- cbind(spa[, 2], spa[, 1])
sc <- myspca2$li   # Only the selected sPCA axes
spe.dh <- dist(sc)
      # Ward's Minimum Variance Clustering
spe.h.ward <- hclust(spe.dh, method = "ward.D")
color <- c("deeppink3", "chartreuse4", "darkgoldenrod2", "darkorange3", "darkorchid4")
# Pour ajouter une forme aux différentes couleurs:
shape <- c(21:25)

k = 5
spebc.ward.g <- cutree(spe.h.ward, k)

# Plot of the Ward clusters on a map of Mikembo
par(mar = c(0, 0, 0, 0))
plot(NULL, xlim = c(-5.1538, -5.1416), ylim = c(33.426, 33.43385))
plot(Azrou_shp, col = "lightgrey", add = T)
grw <- spebc.ward.g
k <- length(levels(factor(grw)))
for(i in 1:k) {
   points(spa[grw == i, 1], spa[grw == i, 2], pch = shape[i], cex = 1.2, col = color[i], 
       bg = color[i])
}
scalebar(0.1, xy = c(-5.143, 33.4265), type = "bar", divs = 2, below = "kilometres")



