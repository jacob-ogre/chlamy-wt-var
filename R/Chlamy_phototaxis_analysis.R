# Calculations for phototaxis quantitative genetics and related.
# Copyright (C) 2013 Kyle Hernandez

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

# Updated July 2013, J. Malcom
# Updated Jan 2014, J. Malcom jmalcom@uconn.edu

library(ggplot2)
library(plyr)
library(lme4)
source("quantgenetics.r")

###########################################################################
# Functions
###########################################################################
subsample <- function(df) {
    subs <- sample(rownames(df), 1)
    newd <- subset(df, rownames(df) != subs)
    return(newd)
}

jackknife <- function(df) {
    Vgs <- rep(0.0, dim(df)[[1]])
    Ves <- rep(0.0, dim(df)[[1]])
    Vps <- rep(0.0, dim(df)[[1]])
    H2s <- rep(0.0, dim(df)[[1]])
    for(i in rownames(df)) {
        newdf <- subset(df, rownames(df) != i)
        res <- quantgene(newdf)
        Vgs[i] <- res[[1]][1]
        Ves[i] <- res[[2]][1]
        Vps[i] <- res[[3]][1]
        H2s[i] <- res[[4]][1]
    }
    boot_res <- data.frame(Vgs, Ves, Vps, H2s)
    return(boot_res)
}

conf_limits <- function(df) {
    for(i in 1:length(df)) {
        print(colnames(df)[i])
        print(quantile(df[,i], probs=c(0.05, 0.5, 0.95)))
    }
}


###########################################################################
# Load data; df name is "d1
###########################################################################
base <- "~/Dropbox/Chlamy_project/Chlamy_phototaxis/"
file <- paste(base, "wt_phototaxis_assays_Nov2012/", 
              "wild_type_phototaxis_filtered.rdata",
              sep="")
load(file)

###########################################################################
# Analysis
###########################################################################
# ANOVA
fit.1 <- lm(d1$INDEX ~ d1$ID)
summary(fit.1)
summary.aov(fit.1)

# Estimate variance components
fit.a <- lmer(d1$INDEX ~ 1 + (1|d1$ID))
summary(fit.a)
Vg <- 0.047845
Ve <- 0.024785
Vp <- Vg + Ve
H <- Vg/Vp
print(c(Vg, Ve, Vp, H))

# Bootstrapping
Boot_res_df <- jackknife(d1)
conf_limits(boot_res_df)

###########################################################################
# Plot the results
###########################################################################
# Estimate means and standard errors of each sample
# Sort based on mean
dd <- ddply(d1, "ID", summarise, mindex = mean(INDEX), sdindex = sd(INDEX))
dd$ID <- reorder(dd$ID, dd$mindex)

# Create dodge and ylim parameters for error bars
dodge <- position_dodge(width=0.9)
limits <- aes(ymax = dd$mindex + dd$sdindex, ymin=dd$mindex - dd$sdindex)

# Plot
pdf(paste, "~/Projects/Chlamy/WT_Natural_Variation/results/figures/",
    "phenotypes/phototaxis_bar.pdf", sep="")

ggplot(dd, aes(x=ID, y=mindex)) +
  geom_bar(stat="identity", position=dodge, col="black", fill="grey60") + 
  geom_errorbar(limits, position=dodge, width=0.25) +
  ylim(-0.5, 1) +
  ylab("Phototaxis Index (SD)") +
  xlab(NULL) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle=90, size=7, vjust=0.5))

dev.off()


