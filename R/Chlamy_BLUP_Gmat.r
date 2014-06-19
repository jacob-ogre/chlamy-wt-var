# Chlamy_wt_fit_var_PCA.r
# PCA decompose the genotypic (breeding) values of Chlamy data.
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

library("FactoMineR")
# setwd(paste("~/Google Drive/Chlamy_project/Chlamy_wt_fitness_assays/", 
#             "extracted_data/", sep=""))

##########################################################################
# Functions
###########################################################################
boot_pca <- function(mat, reps) {
    obs <- scale(mat)
    opc <- PCA(obs, ncp=15, graph=F)
    len_vec <- min(length(obs[1,]), length(obs[,1]))
    refer <- opc$eig$eigenvalue
    for(i in 1:(reps-1)) {
        cur_mat <- apply(mat, 2, sample)
        scaled <- scale(cur_mat)
        pca <- PCA(scaled, ncp=15, graph=F)
        res <- pca$eig$eigenvalue >= opc$eig$eigenvalue
        refer <- cbind(refer, res)
    }
    return(data.frame(refer))
}

cv <- function(x){
    return (100*sd(x) / mean(x))
}

eff_evolv <- function(x) {
    tmp <- sum(x) / max(x)
    return(tmp)
}

max_evolv <- function(x) {
    tmp <- sqrt(max(x))
    return(tmp)
}

###########################################################################
# Load data
###########################################################################
# gcor.dat <- read.table("wt_fitness_maxF_FINAL.tab", sep="\t", header=TRUE)
base <- "/Users/jacobmalcom/Dropbox/Chlamy_project/Chlamy_wt_fitness_assays/extracted_data/"
dat <- read.table(paste(base, "GxE_blups_mod4_mat.tab", sep=""),
                  sep="\t", header=TRUE)
rownames(dat) <- dat[,1]
dat <- dat[,-1]
names(dat)
dat <- dat[,-25]
dat <- dat[,-24]
dat <- dat[,-7]
dat <- dat[,-4]
dat <- dat[,-3]
dim(dat)
names(dat)

means <- read.table(paste(base, "Chlamy_wt_max_linemeans.csv", sep=""),
                      sep=",", header=TRUE)
keep <- c("1373", "1690", "2246", "2290", "2342", "2343", "2344", "2607",
          "2608", "2931", "2932", "2935", "2936", "2937", "2938", "4414",
          "89", "90")
means <- droplevels(subset(means, means$line %in% keep))
rownames(means) <- means[,1]
means <- means[,-1]
names(means)
means <- means[,-25]
means <- means[,-24]
means <- means[,-7]
means <- means[,-4]
means <- means[,-3]
dim(means)

###############################################################################
# plots between line means and blups
# The G-matrix results are very different between using line means and BLUPs,
# so I want to see how the values compare to one-another in graphical form
###############################################################################
par(mfrow=c(5,6))
all_blup <- c()
all_mean <- c()
for(i in 1:length(dat)) {
    plot(abs(dat[,i]) ~ means[,i],
         xlab=paste("Means", names(means)[i], sep=" "),
         ylab=paste("BLUPs", names(dat)[i], sep=" "))
    abline(lm(abs(dat[,1]) ~ means[,i]), col="red")
    all_blup <- c(all_blup, dat[,i])
    all_mean <- c(all_mean, means[,i])
}

par(mfrow=c(1,1))
plot(abs(all_blup) ~ scale(all_mean))
abline(lm(abs(all_blup) ~ scale(all_mean)))

###########################################################################
# PCA of genotypic values
###########################################################################
mean_cent <- scale(t(dat))
cov_dat <- cov(mean_cent)

cov_dat_eig <- eigen(cov_dat)
cov_dat_var <- cov_dat_eig$values / sum(cov_dat_eig$values)
cov_dat_eig_val <- cov_dat_eig$values
cov_dat_eig_val
cov_dat_var

# calculate Kirkpatrick metrics
cov_eig_eff_dims <- eff_evolv(cov_dat_eig$values)
cov_eig_eff_dims

cov_eig_max_evolv <- max_evolv(cov_dat_eig$values)
cov_eig_max_evolv

cov_tot_gen_var <- sum(cov_dat_eig_val)
cov_dat_eig_val[1] * cov_eig_eff_dims
cov_tot_gen_var 

###########################################################################
# bootstrap the eigenvalues
###########################################################################
boots <- boot_pca(mean_cent, 9999)
rowSums(boots) / length(boots)

###########################################################################
# Some plotting
###########################################################################
cov_PCA <- PCA(mean_cent, graph=FALSE)
plot(cov_PCA, choix="ind", new.plot=FALSE)
plot(cov_PCA, choix="var", new.plot=FALSE)

# The FactoMineR biplots can't be combined, at least not in the way I want,
# so will use biplot with a prcomp analysis. The numbers are slightly different
# relative to eigen decomp and FactoMineR (the numbers are the same for these
# two), but because this is just about visualization, and the relative positions
# of the genotypes and environments are effectively identical with the exception
# of a 90 degree rotation, this works.
cov_prc <- prcomp(mean_cent)
biplot(cov_prc, 
       col=c("black", "darkgray"), 
       cex=0.8, 
       xlim=c(-0.5, 0.5),
       ylim=c(-0.5, 0.5),
       xlab="PC1 (27.45%)",
       ylab="PC2 (19.2%)")
