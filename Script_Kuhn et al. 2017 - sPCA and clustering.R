######################################################
#******# Spatial principal component analysis #******#
######################################################

library(adegenet)
library(spdep)
library(adehabitat)
library(SoDA)
library(gclus)
library(cluster)


### I. Queens: L1 and L2 considered together:
#############################################

spe <- read.table("spe_R1_R2_na.txt", h = T, sep = "\t", row.names = 1)
spa1 <- read.table("spa_R1_na.txt", h = T, sep = "\t", row.names = 1)
spa2 <- read.table("spa_R2_na.txt", h = T, sep = "\t", row.names = 1)

SPA <- rbind(spa1, spa2)

spa <- geoXY(SPA[,1], SPA[,2], lat0 = 33.425847, lon0 = -5.154635, unit = 1)

(spegenind <- new(Class = "genind", tab = spe, pop = NULL, ploidy = 2,
   type = "codom"))

spegenind$tab <- scaleGen(spegenind, center = FALSE, scale = FALSE, NA.method = "mean")

myspca <- spca(spegenind, xy = spa, ask = TRUE, plot.nb = TRUE, 
   truenames = TRUE, edit.nb = FALSE)

# Pour génération des figures avec latitudes (ne pas utiliser geoXY) :
# ********************************************************************
# myspca$xy[, 1] <- as.numeric(SPA[, 2])
# myspca$xy[, 2] <- as.numeric(SPA[, 1])

# Proportion de la variabilité décrite par l'axe 1 :

sum <- summary(myspca)$spca[, 2]
# Axe 1 :
sum[1] / sum(sum)

# L'arg. matWeight permet ajouter une matrice de pondération spatiale (= W ?).
# L'arg. scale permet de centrer (ou non) les allèles à variance de 1 (TRUE).
# L'arg. scannf permet de choisir si les valeurs propres doivent être choisies
# interactivement (manuellement) ou pas.
# L'arg. ask permet de décider si les graphes devraient être choisis 
# interactivement ou pas.
# L'arg. plot.nb permet de décider si la figure résultante doit être présentée
# L'arg. truenames (T ou F) --> vrais noms de 'spe' ou labels génériques
# L'arg. edit.nb --> pouvoir éditer la figure résultante manuellement, ou pas.
# L'arg. cn permet d'implémenter notre propre matrice W au lieu du connection
# network.
# Autre arg. -----> ?spca

print(myspca)
summary(myspca, printres = T)

barplot(myspca$eig, main="sPCA eigenvalues", col=spectral(length(myspca$eig)))
legend("topright", fill=spectral(2), leg=c("Global structures", "Local structures"))
abline(h=0,col="grey")

head(myspca$c1)   # Axes de l'analyse stockés en colonnes dans $c1
head(myspca$li)   # Entity scores (coord. des colonies) dans $li

screeplot(myspca)   # Les eigenvalues sont construites à partir de la variance
# ET de l'indice de Moran (I) --> considérer les deux composantes est 
# essentiel pour choisir les axes intéressants à considérer, spatialement. 
# En effet, si grande variance mais faible I, pas de struct. spatiale. Si 
# très faible variance mais I élevé --> probablement même pas de structure 
# génétique.

# Utilisons à présent test de significativité pour décider si les structures
# globales et locales méritent considération.

(myGtest <- global.rtest(spegenind$tab, myspca$lw, nperm=9999))
plot(myGtest)

source("sw.value.R")
sc <- as.vector(myspca$li)[,1]
sw.value(spa, sc)


### II. Queens: L1 alone:
#########################

spe <- read.table("spe_R1_na.txt", h = T, sep = "\t", row.names = 1)
spa <- read.table("spa_R1_na.txt", h = T, sep = "\t", row.names = 1)

spa <- geoXY(spa[,1], spa[,2], lat0 = 33.425847, lon0 = -5.154635, unit = 1)

(spegenind <- new(Class = "genind", tab = spe, pop = NULL, ploidy = 2,
   type = "codom"))

spegenind$tab <- scaleGen(spegenind, center = FALSE, scale = FALSE, NA.method = "mean")

myspca <- spca(spegenind, xy = spa, ask = TRUE, plot.nb = TRUE, 
   truenames = TRUE, edit.nb = FALSE)

# Pour génération des figures avec latitudes (ne pas utiliser geoXY) :
# ********************************************************************
# spa <- read.table("spa_R1_na.txt", h = T, sep = "\t", row.names = 1)
# myspca$xy[, 1] <- as.numeric(spa[, 2])
# myspca$xy[, 2] <- as.numeric(spa[, 1])

# L'arg. matWeight permet ajouter une matrice de pondération spatiale (= W ?).
# L'arg. scale permet de centrer (ou non) les allèles à variance de 1 (TRUE).
# L'arg. scannf permet de choisir si les valeurs propres doivent être choisies
# interactivement (manuellement) ou pas.
# L'arg. ask permet de décider si les graphes devraient être choisis 
# interactivement ou pas.
# L'arg. plot.nb permet de décider si la figure résultante doit être présentée
# L'arg. truenames (T ou F) --> vrais noms de 'spe' ou labels génériques
# L'arg. edit.nb --> pouvoir éditer la figure résultante manuellement, ou pas.
# L'arg. cn permet d'implémenter notre propre matrice W au lieu du connection
# network.
# Autre arg. -----> ?spca

print(myspca)
summary(myspca, printres = T)

# Proportion de la variabilité génétique totale décrite par les trois premiers axes :

sum <- summary(myspca)$spca[, 2]
# Axe 1 :
sum[1] / sum(sum)
# Axe 2 :
sum[2] / sum(sum)
# Axe 3 :
sum[3] / sum(sum)

barplot(myspca$eig, main="sPCA eigenvalues", col=spectral(length(myspca$eig)))
legend("topright", fill=spectral(2), leg=c("Global structures", "Local structures"))
abline(h=0,col="grey")

head(myspca$c1)   # Axes de l'analyse stockés en colonnes dans $c1
head(myspca$li)   # Entity scores (coord. des colonies) dans $li

screeplot(myspca)   # Les eigenvalues sont construites à partir de la variance
# ET de l'indice de Moran (I) --> considérer les deux composantes est 
# essentiel pour choisir les axes intéressants à considérer, spatialement. 
# En effet, si grande variance mais faible I, pas de struct. spatiale. Si 
# très faible variance mais I élevé --> probablement même pas de structure 
# génétique.

# Utilisons à présent test de significativité pour décider si les structures
# globales et locales méritent considération.

(myGtest <- global.rtest(spegenind$tab, myspca$lw, nperm=9999))
plot(myGtest)

source("sw.value.R")
sc <- as.vector(myspca$li)[,1]
sw.value(spa, sc)

#plot(myspca, axis = 1, useLag = T)
colorplot(myspca, axes = 1:3)

(myLtest <- local.rtest(spegenind$tab, myspca$lw, nperm=9999))
plot(myLtest)

   # Clustering of the sPCA axes:
   # ****************************

sc <- myspca$li   # Only the selected sPCA axes

spe.dh <- dist(sc)

      # Ward's Minimum Variance Clustering
spe.h.ward <- hclust(spe.dh, method = "ward.D")

plot(spe.h.ward, cex = 0.6)

         # Looking for interpretable Clusters (Where should the tree be cut?)
         # *********************************************************************        
            # Graph of the Fusion Level Values
            # ********************************

plot(spe.h.ward$height, nrow(sc):2, type = "S", 
   main = "Fusion levels - Hellinger - Ward", ylab = "k (number of clusters)",
   xlab = "h (node height)", col = "grey")
text(spe.h.ward$height, nrow(sc):2, nrow(sc):2, col = "red", cex = 0.8)

            # Graphs of Silhouette Widths
            # ***************************

asw=numeric(nrow(sc))
for(k in 2:(nrow(sc)-1)){
   sil <- silhouette(cutree(spe.h.ward, k = k), spe.dh)
   asw[k] <- summary(sil)$avg.width
}  
k.best <- which.max(asw)   # Best (largest) silhouette width
plot(1: nrow(sc), asw, type = "h" ,main = "Silhouette-optimal number of clusters, 
   Ward", xlab = "k (number of groups)", ylab = "Average silhouette width")
axis(1, k.best, paste("optimum", k.best, sep = "\n"), col = "red", font = 2, 
   col.axis = "red")
points(k.best, max(asw), pch = 16, col = "red", cex = 1.5)

cat("", "Silhouette-optimal number of clusters k = ", k.best, "\n", 
   "with an average silhouette width of", max(asw), "\n")

           # Mantel comparison 
           # *****************
# (Comparison b/ the Distance Matrix and Binary Matrices representing partitions)

# Function to compute a binary distance matrix from groups
grpdist <- function(X)
   {
   require(cluster)
   gr <- as.data.frame(as.factor(X))
   distgr <- daisy(gr, "gower")
   distgr
}
# Run based on the Ward clustering
kt <- data.frame(k = 1: nrow(sc), r = 0)

for(i in 2:(nrow(sc)-1)) {
   gr <- cutree(spe.h.ward, i)
   distgr <- grpdist(gr)
   mt <- cor(spe.dh, distgr, method = "pearson")
   kt[i, 2] <- mt
}
kt
k.best <- which.max(kt$r)
plot(kt$k, kt$r, type = "h", main = "Mantel-optimal number if clusters - Ward", 
   xlab = "k (number of groups)", ylab = "Pearson's correlation")
axis(1, k.best, paste("optimum", k.best, sep="\n"), col = "red", font = 2,
   col.axis = "red")
points(k.best, max(kt$r), pch = 16, col = "red", cex = 1.5)

            # Silhouette Plot of the Final Partition
            # **************************************

# We proceed now to examinate if the group memberships are appropriate
k <- 5   # Nb of clusters
cutg <- cutree(spe.h.ward,k = k)
sil <- silhouette(cutg, spe.dh)
silo <- sortSilhouette(sil)
rownames(silo) <- row.names(sc)[attr(silo, "iOrd")]
plot(silo, main = "Silhouette plot - Hellinger - Ward", cex.names = 0.8,
   col = cutg + 1, nmax.lab = 100)

            # Final Dendrogram with Graphical Options
            # ***************************************

# reorder.hclust reorders objects so that their order in the disimilarity matrix is 
# respected as much as possible.
spe.hwo <- reorder.hclust(spe.h.ward, spe.dh)   
plot(spe.hwo, hang = -1, xlab = "5 groups", sub = "", ylab = "Height", 
   main = "Euclidean - Ward (reordered)", labels = cutree(spe.hwo, k = k), cex = 0.6)
rect.hclust(spe.hwo, k = k)
# Plot the final dendrogram with group colors
source("hcoplot.R")
hcoplot(spe.h.ward, spe.dh, k = 5)   # Adapter la valeur de k ici

            # Spatial Plot of the Clustering Result
            # *************************************

spebc.ward.g <- cutree(spe.h.ward, k)

# Plot of the Ward clusters on a map of Mikembo
plot(spa, asp = 1, type = "n", main = "Five Ward groups", xlab = "x coordinate (m)",
   ylab = "y coordinate (m)")
grw <- spebc.ward.g
k <- length(levels(factor(grw)))
for(i in 1:k) {
   points(spa[grw == i, 1], spa[grw == i, 2], pch = i + 19, cex = 3, col = i + 1,
       bg = i + 1)
}
text(spa, row.names(spa), cex = 0.4, col = "white", font = 2)
legend("bottomright", paste("Group", 1:k), pch = (1:k) + 19, col = 2:(k+1),
   pt.bg = 2:(k+1), pt.cex = 2, bty = "n")

write.table(grw, "grw.QueensL1.txt", sep = "\t")


### III. Queens: L2 alone:
#########################

spe <- read.table("spe_R2_na.txt", h = T, sep = "\t", row.names = 1)
spa <- read.table("spa_R2_na.txt", h = T, sep = "\t", row.names = 1)

spa <- geoXY(spa[,1], spa[,2], lat0 = 33.425847, lon0 = -5.154635, unit = 1)

(spegenind <- new(Class = "genind", tab = spe, pop = NULL, ploidy = 2,
   type = "codom"))

spegenind$tab <- scaleGen(spegenind, center = FALSE, scale = FALSE, NA.method = "mean")

myspca <- spca(spegenind, xy = spa, ask = TRUE, plot.nb = TRUE, 
   truenames = TRUE, edit.nb = FALSE)

# Pour génération des figures avec latitudes (ne pas utiliser geoXY) :
# ********************************************************************
# spa <- read.table("spa_R2_na.txt", h = T, sep = "\t", row.names = 1)
# myspca$xy[, 1] <- as.numeric(spa[, 2])
# myspca$xy[, 2] <- as.numeric(spa[, 1])

# L'arg. matWeight permet ajouter une matrice de pondération spatiale (= W ?).
# L'arg. scale permet de centrer (ou non) les allèles à variance de 1 (TRUE).
# L'arg. scannf permet de choisir si les valeurs propres doivent être choisies
# interactivement (manuellement) ou pas.
# L'arg. ask permet de décider si les graphes devraient être choisis 
# interactivement ou pas.
# L'arg. plot.nb permet de décider si la figure résultante doit être présentée
# L'arg. truenames (T ou F) --> vrais noms de 'spe' ou labels génériques
# L'arg. edit.nb --> pouvoir éditer la figure résultante manuellement, ou pas.
# L'arg. cn permet d'implémenter notre propre matrice W au lieu du connection
# network.
# Autre arg. -----> ?spca

print(myspca)
summary(myspca, printres = T)

# Proportion de la variabilité décrite par les trois premiers axes :

sum <- summary(myspca)$spca[, 2]
# Axe 1 :
sum[1] / sum(sum)
# Axe 2 :
sum[2] / sum(sum)
# Axe 3 :
sum[3] / sum(sum)

barplot(myspca$eig, main="sPCA eigenvalues", col=spectral(length(myspca$eig)))
legend("topright", fill=spectral(2), leg=c("Global structures", "Local structures"))
abline(h=0,col="grey")

head(myspca$c1)   # Axes de l'analyse stockés en colonnes dans $c1
head(myspca$li)   # Entity scores (coord. des colonies) dans $li

screeplot(myspca)   # Les eigenvalues sont construites à partir de la variance
# ET de l'indice de Moran (I) --> considérer les deux composantes est 
# essentiel pour choisir les axes intéressants à considérer, spatialement. 
# En effet, si grande variance mais faible I, pas de struct. spatiale. Si 
# très faible variance mais I élevé --> probablement même pas de structure 
# génétique.

# Utilisons à présent test de significativité pour décider si les structures
# globales et locales méritent considération.

(myGtest <- global.rtest(spegenind$tab, myspca$lw, nperm=9999))
plot(myGtest)

source("sw.value.R")
sc <- as.vector(myspca$li)[,1]
sw.value(spa, sc)

#plot(myspca, axis = 1, useLag = T)
colorplot(myspca, axes = 1:3)

(myLtest <- local.rtest(spegenind$tab, myspca$lw, nperm=9999))
plot(myLtest)

   # Clustering of the sPCA axes:
   # ****************************

sc <- myspca$li   # Only the selected sPCA axes

spe.dh <- dist(sc)

      # Ward's Minimum Variance Clustering
spe.h.ward <- hclust(spe.dh, method = "ward.D")

plot(spe.h.ward, cex = 0.6)

         # Looking for interpretable Clusters (Where should the tree be cut?)
         # *********************************************************************        
            # Graph of the Fusion Level Values
            # ********************************

plot(spe.h.ward$height, nrow(sc):2, type = "S", 
   main = "Fusion levels - Hellinger - Ward", ylab = "k (number of clusters)",
   xlab = "h (node height)", col = "grey")
text(spe.h.ward$height, nrow(sc):2, nrow(sc):2, col = "red", cex = 0.8)

            # Graphs of Silhouette Widths
            # ***************************

asw=numeric(nrow(sc))
for(k in 2:(nrow(sc)-1)){
   sil <- silhouette(cutree(spe.h.ward, k = k), spe.dh)
   asw[k] <- summary(sil)$avg.width
}  
k.best <- which.max(asw)   # Best (largest) silhouette width
plot(1: nrow(sc), asw, type = "h" ,main = "Silhouette-optimal number of clusters, 
   Ward", xlab = "k (number of groups)", ylab = "Average silhouette width")
axis(1, k.best, paste("optimum", k.best, sep = "\n"), col = "red", font = 2, 
   col.axis = "red")
points(k.best, max(asw), pch = 16, col = "red", cex = 1.5)

cat("", "Silhouette-optimal number of clusters k = ", k.best, "\n", 
   "with an average silhouette width of", max(asw), "\n")

           # Mantel comparison 
           # *****************
# (Comparison b/ the Distance Matrix and Binary Matrices representing partitions)

# Function to compute a binary distance matrix from groups
grpdist <- function(X)
   {
   require(cluster)
   gr <- as.data.frame(as.factor(X))
   distgr <- daisy(gr, "gower")
   distgr
}
# Run based on the Ward clustering
kt <- data.frame(k = 1: nrow(sc), r = 0)

for(i in 2:(nrow(sc)-1)) {
   gr <- cutree(spe.h.ward, i)
   distgr <- grpdist(gr)
   mt <- cor(spe.dh, distgr, method = "pearson")
   kt[i, 2] <- mt
}
kt
k.best <- which.max(kt$r)
plot(kt$k, kt$r, type = "h", main = "Mantel-optimal number if clusters - Ward", 
   xlab = "k (number of groups)", ylab = "Pearson's correlation")
axis(1, k.best, paste("optimum", k.best, sep="\n"), col = "red", font = 2,
   col.axis = "red")
points(k.best, max(kt$r), pch = 16, col = "red", cex = 1.5)

            # Silhouette Plot of the Final Partition
            # **************************************

# We proceed now to examinate if the group memberships are appropriate
k <- 5   # Nb of clusters
cutg <- cutree(spe.h.ward,k = k)
sil <- silhouette(cutg, spe.dh)
silo <- sortSilhouette(sil)
rownames(silo) <- row.names(sc)[attr(silo, "iOrd")]
plot(silo, main = "Silhouette plot - Hellinger - Ward", cex.names = 0.8,
   col = cutg + 1, nmax.lab = 100)

            # Final Dendrogram with Graphical Options
            # ***************************************

# reorder.hclust reorders objects so that their order in the disimilarity matrix is 
# respected as much as possible.
spe.hwo <- reorder.hclust(spe.h.ward, spe.dh)   
plot(spe.hwo, hang = -1, xlab = "5 groups", sub = "", ylab = "Height", 
   main = "Euclidean - Ward (reordered)", labels = cutree(spe.hwo, k = k), cex = 0.6)
rect.hclust(spe.hwo, k = k)
# Plot the final dendrogram with group colors
source("hcoplot.R")
hcoplot(spe.h.ward, spe.dh, k = 5)   # Adapter la valeur de k ici

            # Spatial Plot of the Clustering Result
            # *************************************

spebc.ward.g <- cutree(spe.h.ward, k)

# Plot of the Ward clusters on a map of Mikembo
plot(spa, asp = 1, type = "n", main = "Five Ward groups", xlab = "x coordinate (m)",
   ylab = "y coordinate (m)")
grw <- spebc.ward.g
k <- length(levels(factor(grw)))
for(i in 1:k) {
   points(spa[grw == i, 1], spa[grw == i, 2], pch = i + 19, cex = 3, col = i + 1,
       bg = i + 1)
}
text(spa, row.names(spa), cex = 0.4, col = "white", font = 2)
legend("bottomright", paste("Group", 1:k), pch = (1:k) + 19, col = 2:(k+1),
   pt.bg = 2:(k+1), pt.cex = 2, bty = "n")

write.table(grw, "grw.QueensL2.txt", sep = "\t")


### IV. Males: SP1 alone:
#########################

spe <- read.table("spe_SP1_na.txt", h = T, sep = "\t", row.names = 1)
spa <- read.table("spa_SP1_na.txt", h = T, sep = "\t", row.names = 1)

spa <- geoXY(spa[,1], spa[,2], lat0 = 33.425847, lon0 = -5.154635, unit = 1)

spax <- jitter(spa[,1], factor=2)
spay <- jitter(spa[,2], factor=2)
spa <- cbind(spax, spay)

erase <- as.vector(which(colSums(spe)==0))   # Enlever allèles inutiles
spe <- spe[,-erase]

(spegenind <- new(Class = "genind", tab = spe, pop = NULL, ploidy = 2,
   type = "codom"))

spegenind$tab <- scaleGen(spegenind, center = FALSE, scale = FALSE, NA.method = "mean")

myspca <- spca(spegenind, xy = spa, ask = TRUE, plot.nb = TRUE, 
   truenames = TRUE, edit.nb = FALSE)

# Pour génération des figures avec latitudes (ne pas utiliser geoXY) :
# ********************************************************************
# spa <- read.table("spa_SP1_na.txt", h = T, sep = "\t", row.names = 1)
# myspca$xy[, 1] <- as.numeric(spa[, 2])
# myspca$xy[, 2] <- as.numeric(spa[, 1])

# L'arg. matWeight permet ajouter une matrice de pondération spatiale (= W ?).
# L'arg. scale permet de centrer (ou non) les allèles à variance de 1 (TRUE).
# L'arg. scannf permet de choisir si les valeurs propres doivent être choisies
# interactivement (manuellement) ou pas.
# L'arg. ask permet de décider si les graphes devraient être choisis 
# interactivement ou pas.
# L'arg. plot.nb permet de décider si la figure résultante doit être présentée
# L'arg. truenames (T ou F) --> vrais noms de 'spe' ou labels génériques
# L'arg. edit.nb --> pouvoir éditer la figure résultante manuellement, ou pas.
# L'arg. cn permet d'implémenter notre propre matrice W au lieu du connection
# network.
# Autre arg. -----> ?spca

print(myspca)
summary(myspca, printres = T)

# Proportion de la variabilité décrite par les trois premiers axes :

sum <- summary(myspca)$spca[, 2]
# Axe 1 :
sum[1] / sum(sum)
# Axe 2 :
sum[2] / sum(sum)
# Axe 3 :
sum[3] / sum(sum)

barplot(myspca$eig, main="sPCA eigenvalues", col=spectral(length(myspca$eig)))
legend("topright", fill=spectral(2), leg=c("Global structures", "Local structures"))
abline(h=0,col="grey")

head(myspca$c1)   # Axes de l'analyse stockés en colonnes dans $c1
head(myspca$li)   # Entity scores (coord. des colonies) dans $li

screeplot(myspca)   # Les eigenvalues sont construites à partir de la variance
# ET de l'indice de Moran (I) --> considérer les deux composantes est 
# essentiel pour choisir les axes intéressants à considérer, spatialement. 
# En effet, si grande variance mais faible I, pas de struct. spatiale. Si 
# très faible variance mais I élevé --> probablement même pas de structure 
# génétique.

# Utilisons à présent test de significativité pour décider si les structures
# globales et locales méritent considération.

(myGtest <- global.rtest(spegenind$tab, myspca$lw, nperm=999))
plot(myGtest)

source("sw.value.R")
sc <- as.vector(myspca$li)[,1]
sw.value(spa, sc)

#plot(myspca, axis = 1, useLag = T)
colorplot(myspca, axes = 1:3)

(myLtest <- local.rtest(spegenind$tab, myspca$lw, nperm=999))
plot(myLtest)

### V. Males: SP2 alone:
########################

spe <- read.table("spe_SP2_na.txt", h = T, sep = "\t", row.names = 1)
spa <- read.table("spa_SP2_na.txt", h = T, sep = "\t", row.names = 1)

spa <- geoXY(spa[,1], spa[,2], lat0 = 33.425847, lon0 = -5.154635, unit = 1)

spax <- jitter(spa[,1], factor=2)
spay <- jitter(spa[,2], factor=2)
spa <- cbind(spax, spay)

(spegenind <- new(Class = "genind", tab = spe, pop = NULL, ploidy = 2,
   type = "codom"))

spegenind$tab <- scaleGen(spegenind, center = FALSE, scale = FALSE, NA.method = "mean")

myspca <- spca(spegenind, xy = spa, ask = TRUE, plot.nb = TRUE, 
   truenames = TRUE, edit.nb = FALSE)

# Pour génération des figures avec latitudes (ne pas utiliser geoXY) :
# ********************************************************************
# spa <- read.table("spa_SP2_na.txt", h = T, sep = "\t", row.names = 1)
# myspca$xy[, 1] <- as.numeric(spa[, 2])
# myspca$xy[, 2] <- as.numeric(spa[, 1])

# L'arg. matWeight permet ajouter une matrice de pondération spatiale (= W ?).
# L'arg. scale permet de centrer (ou non) les allèles à variance de 1 (TRUE).
# L'arg. scannf permet de choisir si les valeurs propres doivent être choisies
# interactivement (manuellement) ou pas.
# L'arg. ask permet de décider si les graphes devraient être choisis 
# interactivement ou pas.
# L'arg. plot.nb permet de décider si la figure résultante doit être présentée
# L'arg. truenames (T ou F) --> vrais noms de 'spe' ou labels génériques
# L'arg. edit.nb --> pouvoir éditer la figure résultante manuellement, ou pas.
# L'arg. cn permet d'implémenter notre propre matrice W au lieu du connection
# network.
# Autre arg. -----> ?spca

print(myspca)
summary(myspca, printres = T)

# Proportion de la variabilité décrite par les trois premiers axes :

sum <- summary(myspca)$spca[, 2]
# Axe 1 :
sum[1] / sum(sum)
# Axe 2 :
sum[2] / sum(sum)
# Axe 3 :
sum[3] / sum(sum)

barplot(myspca$eig, main="sPCA eigenvalues", col=spectral(length(myspca$eig)))
legend("topright", fill=spectral(2), leg=c("Global structures", "Local structures"))
abline(h=0,col="grey")

head(myspca$c1)   # Axes de l'analyse stockés en colonnes dans $c1
head(myspca$li)   # Entity scores (coord. des colonies) dans $li

screeplot(myspca)   # Les eigenvalues sont construites à partir de la variance
# ET de l'indice de Moran (I) --> considérer les deux composantes est 
# essentiel pour choisir les axes intéressants à considérer, spatialement. 
# En effet, si grande variance mais faible I, pas de struct. spatiale. Si 
# très faible variance mais I élevé --> probablement même pas de structure 
# génétique.

# Utilisons à présent test de significativité pour décider si les structures
# globales et locales méritent considération.

(myGtest <- global.rtest(spegenind$tab, myspca$lw, nperm=999))
plot(myGtest)

source("sw.value.R")
sc <- as.vector(myspca$li)[,1]
sw.value(spa, sc)

plot(myspca, axis = 1, useLag = T)
colorplot(myspca, axes = 1:3)

(myLtest <- local.rtest(spegenind$tab, myspca$lw, nperm=999))
plot(myLtest)
