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

###########################################################################
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
base <- "/home/jacob/Dropbox/Chlamy_project/Chlamy_wt_fitness_assays/extracted_data/"
dat <- read.table(paste(base, "GxE_blups_mod4_mat.tab", sep=""),
                  sep="\t", header=TRUE)
head(dat)

# some df management, including remove N++100, which is synonymous with N++
rownames(dat) <- dat[,1]
dat <- dat[,-1]
dat <- dat[,-7]

###########################################################################
# PCA of genotypic values
###########################################################################
mean_cent <- (t(dat) - colMeans(t(dat), na.rm=T))
mean_cent <- mean_cent  / apply(mean_cent, 2, sd, na.rm=T)
cov_dat <- cov(mean_cent, use="pairwise")
cor_dat <- cor(mean_cent, use="pairwise")

cov_dat_eig <- eigen(cov_dat)
cov_dat_var <- cov_dat_eig$values / sum(cov_dat_eig$values)
cor_dat_eig <- eigen(cor_dat)
cor_dat_var <- cor_dat_eig$values / sum(cor_dat_eig$values)
cov_dat_var
cor_dat_var

# calculated eff
cov_eig_eff_dims <- eff_evolv(cov_dat_eig$values)
cor_eig_eff_dims <- eff_evolv(cor_dat_eig$values)
cov_eig_eff_dims
cor_eig_eff_dims

cov_eig_max_evolv <- max_evolv(cov_dat_eig$values)
cor_eig_max_evolv <- max_evolv(cor_dat_eig$values)
cov_eig_max_evolv
cor_eig_max_evolv

###########################################################################
# clustering
###########################################################################
datamat <- scale(dat[,2:length(dat)])
rownames(datamat) <- gcor.dat[,1]
dist_mat <- dist(datamat, "euclidian")
fit <- hclust(dist_mat, method="ward")

plot(fit, main="", xlab="", ylab="distance")

###########################################################################
# Write results
###########################################################################
write.table(line_fold_df,
            file="wt_fit_var_linewise_fold_differences.tab",
            sep="\t",
            eol="\n",
            row.names=FALSE,
            quote=FALSE)

write.table(env_fold_df,
            file="wt_fit_var_envwise_fold_differences.tab",
            sep="\t",
            eol="\n",
            row.names=FALSE,
            quote=FALSE)

