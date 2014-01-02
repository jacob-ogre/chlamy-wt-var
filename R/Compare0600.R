# Compare size at the 06:00h timepoint.
# Copyright (C) 2014 Jacob Malcom, jmalcom@uconn.edu

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
source("myPlotMeans.topleft.R")

dat <- read.table("Chlamy_cell_sizes_0600comparison.csv", sep=",", header=TRUE)
names(dat)

# convert to factors:
dat$Line <- as.factor(dat$Line)
dat$Time <- as.factor(dat$Time)
dat$Image <- as.factor(dat$Image)

# plot means +- 2se
new_ord <- reorder(dat$Line, dat$Mean, FUN=mean)
myPlotMeans.topleft(dat$Mean, new_ord, dat$Time, xlab="Line", 
                    ylab="Cell Area (um^2)", legend.lab="Time", cex=1)

new_ord <- reorder(dat$Line, dat$Cells, FUN=mean)
myPlotMeans.topleft(dat$Cells, new_ord, dat$Time, xlab="Line", 
                    ylab="Cell Count", legend.lab="Time", cex=1)

new_ord <- reorder(dat$Line, dat$Area, FUN=mean)
myPlotMeans.topleft(dat$Area, new_ord, dat$Time, xlab="Line", 
                    ylab="Cluster Area (um^2)", legend.lab="Time", cex=1)
