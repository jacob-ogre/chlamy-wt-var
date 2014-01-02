# Multipanel of Chlamy sizes.
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

# Create figure of growth comparisons between the 18 lines at 2 time points 06:00 and 12:00 and
# show the growth patterns of the 2 lines with all time points.

# Load sources
library(sciplot)
source("~/Projects/Chlamy/WT_Natural_Variation/bin/chlamy-mapping-population/rscripts/general_plot_fx/myPlotMedian.general.R")

# Read in the 1690 v 2290 all time points file
allTP <- read.table("~/Dropbox/Chlamy_project/Chlamy_photos/Chlamy_1690_2290_comparison2.csv", sep=",", header = TRUE)

# convert to factors:
allTP$Line <- as.factor(allTP$Line)
allTP$Time <- as.factor(allTP$Time)
allTP$Image <- as.factor(allTP$Image)
# Reorder
new_ord <- reorder(allTP$Time, allTP$Order, mean)

# Read in the 18-line 2 timepoint file
twoTP <- read.table("~/Dropbox/Chlamy_project/Chlamy_photos/Chlamy_2_time_point_all.csv", sep=",", header = TRUE)
# convert to factors:
twoTP$Line <- as.factor(twoTP$Line)
twoTP$Time <- as.factor(twoTP$Time)
twoTP$Image <- as.factor(twoTP$Image)

def.par <- par(no.readonly=TRUE)
# Output to pdf
pdf("~/Dropbox/Chlamydomonas/Wild-type_variation/Figures/base_figures/cell-area-count-figure.pdf", 10, 6)

# Make plot area
layout(matrix(c(1,2,3,4,5,6), 2, 3,  byrow=TRUE), widths=c(4,3,3), heights=c(1,1))
# Cell area all TP
par(mar=c(3,3,1,1))
myPlotMedian.general(allTP$Mean, new_ord, allTP$Line, xlab="Time", 
                     ylab=expression(paste(bold("Cell Area "), "(", mu, m^2, ")")), 
                     legend.lab="Line", ylim=c(0,5),
                     cex=0.75, pos="topleft", cex.lab=1.2, error.width=0.05)

# twoTP area bar
par(mar=c(5,5,1,1))
tp1 <- twoTP[which(twoTP$Time=="06:00"),]
bargraph.CI(tp1$Line, tp1$Mean, xlab="",
            cex.names=0.75, las=2, 
            ylab=expression(paste(bold("Cell Area "), "(", mu, m^2, ")")))

# twoTP count bar
par(mar=c(5,5,1,1))
bargraph.CI(tp1$Line, tp1$Cells, xlab="",
            cex.names=0.75, las=2, 
            ylab="Cell Count")

# Side text
mtext("06:00hr", side=4, outer=FALSE, line=0.1, cex = 0.8)

# Cell count all TP
par(mar=c(3,3,1,1))
myPlotMedian.general(allTP$Cells, new_ord, allTP$Line, xlab="Time", 
                     ylab="Cell Count", legend.lab="Line", ylim=c(0,5),
                     cex=0.75, pos="topright", cex.lab = 1.2, error.width=0.05)

# twoTP area bar
par(mar=c(5,5,1,1))
tp2 <- twoTP[which(twoTP$Time=="12:00"),]
bargraph.CI(tp2$Line, tp2$Mean, xlab="",
            cex.names=0.75, las=2, cex.lab=1.2,
            ylab=expression(paste(bold("Cell Area "), "(", mu, m^2, ")")))

# twoTP count bar
par(mar=c(5,5,1,1))
bargraph.CI(tp2$Line, tp2$Cells, xlab="",
            cex.names=0.75, las=2, ylab=paste("Cell Count"))

# Side text
mtext("12:00hr", side=4, outer=FALSE, line=0.1, cex = 0.8)

par(def.par)

dev.off()
