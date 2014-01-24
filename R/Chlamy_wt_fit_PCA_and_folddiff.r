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
setwd(paste("~/Dropbox/Chlamy_project/Chlamy_wt_fitness_assays/", 
            "extracted_data/", sep=""))

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

###########################################################################
# Load data
###########################################################################
# gcor.dat <- read.table("wt_fitness_maxF_FINAL.tab", sep="\t", header=TRUE)
dat <- read.table("wt_fitness_maxF_FINAL.tab", sep="\t", header=TRUE)

# remove N++100, which is synonymous with N++
gcor.dat <- gcor.dat[,-7]

###############################################################################
# calculate line means
###############################################################################
means <- data.frame(t(tapply(dat$Fluor, 
                             INDEX=list(dat$Env, dat$Line), 
                             FUN=mean)))
gcor.dat <- means[,-15]
head(means)

write.table(gcor.dat, 
            file="wt_fitness_max_MEANS_FINAL.tab",
            sep="\t",
            quote=FALSE)

###########################################################################
# PCA of genotypic values
###########################################################################
gen_pca <- PCA(gcor.dat)
dim_des <- dimdesc(gen_pca, axes=1:5)
pca_axs <- data.frame(gen_pca$ind$coord[,1:2])

write.table(pca_axs,
            file="wt_var_Gmat_PCs_FINAL.tab",
            sep="\t",
            quote=FALSE,
            eol="\n")

boot_res <- boot_pca(gcor.dat[,2:length(gcor.dat)], 9999)
rowSums(boot_res[,2:length(boot_res)])

#############################    
# check factoMineR against manual:
a <- cor(means, use="pairwise")
a_eig <- eigen(a)
a_var <- a_eig$values / sum(a_eig$values)

c <- scale(gcor.dat)
d <- cov(c, use="pairwise")
d_eig <- eigen(d)
d_var <- d_eig$values / sum(d_eig$values)

a_var
d_var
gen_pca$eig$`percentage of variance`

###############################################################################
# Plot the PCA
###############################################################################
fig_base <- "~/Dropbox/writings/mss/Chlamydomonas/Wild-type_variation/Figures/"
fig_file <- paste(fig_base, "fitness_PCA.pdf", sep="")

pdf(fig_file, height=10, width=10)

plot.PCA(gen_pca, choix="var")

dev.off()

###########################################################################
# clustering
###########################################################################
datamat <- scale(gcor.dat[,2:length(gcor.dat)])
rownames(datamat) <- gcor.dat[,1]
dist_mat <- dist(datamat, "euclidian")
fit <- hclust(dist_mat, method="ward")

plot(fit, main="", xlab="", ylab="distance")

###########################################################################
# Get fold-change differences since the data is loaded
###########################################################################
line_maxes <- as.numeric(apply(gcor.dat[,2:length(gcor.dat)], MARGIN=1, 
                               FUN=max))
line_mins <- as.numeric(apply(gcor.dat[,2:length(gcor.dat)], MARGIN=1, 
                              FUN=min))
line_cvs  <- as.numeric(apply(gcor.dat[,2:length(gcor.dat)], MARGIN=1, FUN=cv))
line_fold <- line_maxes / line_mins

env_maxes <- as.numeric(apply(gcor.dat[,2:length(gcor.dat)], MARGIN=2, 
                              FUN=max))
env_mins <- as.numeric(apply(gcor.dat[,2:length(gcor.dat)], MARGIN=2, 
                             FUN=min))
env_cvs  <- as.numeric(apply(gcor.dat[,2:length(gcor.dat)], MARGIN=2, FUN=cv))
env_fold <- env_maxes / env_mins

line_fold_df <- data.frame(gcor.dat[,1], line_maxes, line_mins, line_fold,
                           line_cvs)
colnames(line_fold_df) <- c("line", "max", "min", "fold_difference", "CV")
env_fold_df <- data.frame(colnames(gcor.dat)[2:length(gcor.dat)], env_maxes, 
                          env_mins, env_fold, env_cvs)
colnames(env_fold_df) <- c("env",  "max", "min", "fold_difference", "CV")

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

