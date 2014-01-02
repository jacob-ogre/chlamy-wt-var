# Chlamy_1690_2290_size_comparisons.R
# Comparison of cell area, size.

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


library(sciplot)
source("myPlotMeans.general.R")
source("myPlotMedian.general.R")

###########################################################################
# load data
###########################################################################
dat <- read.table("Chlamy_1690_2290_comparison2.csv", sep=",", header=TRUE)
names(dat)

# convert to factors:
dat$Line <- as.factor(dat$Line)
dat$Time <- as.factor(dat$Time)
dat$Image <- as.factor(dat$Image)

# look at ditributions:
hist(dat$Cells)
hist(dat$Mean)
hist(dat$Area)

###########################################################################
# Plot the data
###########################################################################
par(mfrow = c(2,1))
new_ord <- reorder(dat$Time, dat$Order, mean)
myPlotMeans.general(dat$Mean, new_ord, dat$Line, xlab="", 
                    ylab="Cell Area (um^2)", legend.lab="Line", cex=1)

myPlotMedian.general(dat$Cells, new_ord, dat$Line, xlab="Time", 
                    ylab="Cell Count (median)", legend.lab="Line", 
                    ylim=c(0,5), cex=1, pos="topright")

myPlotMeans.general(dat$Area, new_ord, dat$Line, xlab="Time", 
                     ylab="Clump Area (um^2)", legend.lab="Line", 
                     ylim=c(0,5), cex=1, pos="topleft")
