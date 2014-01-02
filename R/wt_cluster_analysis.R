# wt_cluster_analysis.r
# Perform basic euclidian clustering on phenotypic and genotypic data.
# Copyright (C) 2013 Jacob Malcom, jmalcom@uconn.edu

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>.



setwd("~/Dropbox/Chlamy_project/")
phen_bas <- "Chlamy_wt_fitness_assays/wt_fitness_assays_Oct2012/extracted_data/"
geno_bas <- "Chlamy_genomics/"

###########################################################################
# Functions
###########################################################################
normalize <- function (df) {
    norm1 <- df / rowMeans(df, na.rm=TRUE)
    Rsums <- (rowSums(df, na.rm=TRUE) + 1) / (2 + 2 * length(df[1,]))
    norm2 <- norm1 / sqrt(Rsums*(1-Rsums))
    return(norm2)
}

###########################################################################
# Load data
###########################################################################
phen_fil <- paste(phen_bas, "Chlamy_wt_max_linemeans.csv", sep="")
phen_dat <- read.table(phen_fil, sep=",", header=TRUE)
# remove N++100, which is synonymous with N++
phen_dat <- phen_dat[,-7]
rownames(phen_dat) <- phen_dat[,1]

geno_fil <- paste(geno_bas, "Cre_poly75_numeric.tab", sep="")
geno_dat <- read.table(geno_fil, sep="\t", header=TRUE)
rownames(geno_dat) <- geno_dat[,1]
norm_gen <- normalize(geno_dat[,2:length(geno_dat)])

###########################################################################
# clustering
###########################################################################
phen_mat <- scale(phen_dat[,2:length(phen_dat)])
phen_dist <- dist(phen_mat, "euclidian")
phen_fit <- hclust(phen_dist, method="ward")
phen_nj <- nj(phen_dist)

phen_dat2 <- phen_dat[c(-24, -23, -22, -5, -3, -1),]
phen_mat2 <- scale(phen_dat2[,2:length(phen_dat2)])
phen_dist2 <- dist(phen_mat2, "euclidian")
phen_fit2 <- hclust(phen_dist2, method="ward")
phen2_nj <- nj(phen_dist2)

geno_mat <- scale(t(norm_gen))
geno_dist <- dist(geno_mat)
geno_fit <- hclust(geno_dist, method="ward")
geno_nj <- nj(geno_dist)

###########################################################################
# Plots
###########################################################################
par(mfrow=c(1,3))
plot(phen_fit, main="Phenotype Clustering", xlab="", ylab="distance")
plot(phen_fit2, main="Phenotype Clustering", xlab="", ylab="distance")
plot(geno_fit, main="Genotype Clustering", xlab="", ylab="distance")

# going with the nj trees, at least for now
outbase <- "~/Dropbox/writings/mss/Chlamydomonas/Wild-type_variation/Figures/"
outfil <- paste(outbase, "phen_geno_clusters_nj.pdf", sep="")

pdf(outfil, height=6, width=12)

par(mfrow=c(1,3))
plot(phen_nj, font=1, cex=1.5, main="All phenotypes")
plot(phen2_nj, font=1, cex=1.5, main="wt-only phenotypes")
plot(geno_nj, font=1, cex=1.5, main="wt genotypes")

dev.off()
