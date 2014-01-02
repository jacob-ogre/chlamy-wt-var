# Create reaction norm of 18 lines across N, S, and P gradients.
# Copyright (C) 2013 J. Malcom and K. Hernandez

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

# Load sources
source("myPlotMean.general.R")

# Read in MAX file
base <- "~/Dropbox/Chlamy_project/Chlamy_wt_fitness_assays/"
maxf <- paste(base, "extracted_data/Chlamy_wt_assay_max_ODs.csv",
              sep="")
max.2 <- read.table(maxf, sep=",", header=TRUE)

max.2$treatment <- as.factor(max.2$treatment)
max.2$line <- as.factor(max.2$line)

head(max.2)
levels(max.2$treatment)

###########################################################################
# Plot
###########################################################################
outbase <- "~/Dropbox/writings/mss/Chlamydomonas/Wild-type_variation/"
pdf(paste(outbase, "Figures/base_figures/reaction-norm-figure.pdf",
          sep=""), width=12, height=4)

par(mfrow=c(1,3))
####---N gradient!---####
N.grad.dat <- data.frame()
N10 <- subset.data.frame(max.2, max.2$treatment == "N10")
N50 <- subset.data.frame(max.2, max.2$treatment == "N50")
N100 <- subset.data.frame(max.2, max.2$treatment == "N++")
N150 <- subset.data.frame(max.2, max.2$treatment == "N150")

N.grad.dat <- rbind(N10, N50, N100, N150)
N.grad.dat <- droplevels(N.grad.dat)
N.grad.dat$treatment = factor(N.grad.dat$treatment,
                              levels(N.grad.dat$treatment)[c(2,4,1,3)], ordered=T)
# rm(N10, N50, N100, N150)

myPlotMean.general(N.grad.dat$max.measure, N.grad.dat$treatment, 
                   N.grad.dat$line, xlab="Treatment", ylab="OD 440/680 (max)",
                   pos="", legend.lab="Line", cex=1)

####---P gradient!---####
P10 <- subset.data.frame(max.2, max.2$treatment == "P10")
P50 <- subset.data.frame(max.2, max.2$treatment == "P50")
P100 <- subset.data.frame(max.2, max.2$treatment == "N++")
P150 <- subset.data.frame(max.2, max.2$treatment == "P150")
P.grad.dat <- rbind(P10, P50, P100, P150)
P.grad.dat <- droplevels(P.grad.dat)
levels(P.grad.dat$treatment)
P.grad.dat$treatment = factor(P.grad.dat$treatment,
                              levels(P.grad.dat$treatment)[c(2,4,1,3)], ordered=T)
# rm(P10, P50, P150)

myPlotMean.general(P.grad.dat$max.measure, P.grad.dat$treatment, 
                   P.grad.dat$line, xlab="Treatment", ylab="OD 440/680 (max)",
                   legend.lab="line", cex=1, pos="")

####---S gradient!---####
S10 <- subset.data.frame(max.2, max.2$treatment == "S10")
S50 <- subset.data.frame(max.2, max.2$treatment == "S50")
S100 <- subset.data.frame(max.2, max.2$treatment == "N++")
S150 <- subset.data.frame(max.2, max.2$treatment == "S150")
S.grad.dat <- rbind(S10, S50, S100, S150)
S.grad.dat <- droplevels(S.grad.dat)
levels(S.grad.dat$treatment)
S.grad.dat$treatment = factor(S.grad.dat$treatment,
                              levels(S.grad.dat$treatment)[c(2,4,1,3)], ordered=T)
# rm(S10, S50, S150)

myPlotMean.general(S.grad.dat$max.measure, S.grad.dat$treatment, 
                   S.grad.dat$line, xlab="Treatment", ylab="OD 440/680 (max)",
                   legend.lab="line", cex=1, pos="")
dev.off()

