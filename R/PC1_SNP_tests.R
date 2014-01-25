# Test associations between first eigenvector and SNPs.
# Copyright (C) 2014 Jacob Malcom, jmalcom@uconn.edu
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>.

###############################################################################
# Import data
###############################################################################
base <- "~/Dropbox/Chlamy_project/"
PCsf <- paste(base, 
              "Chlamy_wt_fitness_assays/extracted_data/wt_var_Gmat_PCs_FINAL.tab",
              sep="")
SNPf <- paste(base, 
              "Chlamy_wt_geno_mats/Cre75_trans.tab", sep="")

PCs_dat <- read.table(PCsf, sep="\t", header=TRUE)
PCs_dat <- PCs_dat[-16:-21,]
PCs_dat <- PCs_dat[-19,]
rownames(PCs_dat)

SNP_dat <- read.table(SNPf, sep="\t", header=TRUE) 
rownames(SNP_dat) <- SNP_dat[,1]
rownames(SNP_dat)

###############################################################################
# Basic analysis
###############################################################################
res <- rep(0.0, length(SNP_dat) - 1)
for(i in 2:length(SNP_dat)) {
    a_mod <- lm(PCs_dat$Dim.1 ~ SNP_dat[,i-1])
    a_aov <- aov(a_mod)
    res[i-1] <- summary(a_aov)[[1]][5][[1]][1]
}
results <- data.frame(colnames(SNP_dat)[2:length(SNP_dat)], res)
colnames(results) <- c("SNP", "p.val")
nom_sig <- subset(results, results$p.val < 0.05)
nom_sig

outf <- paste(base, "Chlamy_wt_fitness_assays/sig_SNPs.tab", sep="")
write.table(nom_sig,
            outf,
            sep="\t",
            quote=FALSE)
