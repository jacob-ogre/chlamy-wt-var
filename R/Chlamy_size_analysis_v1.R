# Analyze Chlamy wt size variation.
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
source("myPlotMeans.topleft.R")

###########################################################################
# Load data
###########################################################################
## GOING TO HAVE TO CHANGE/SET THE PATH
data <- read.table("Chlamy_cell_sizes_JWMcopy.csv", sep=",", header = TRUE)
names(data)

# convert to factors:
data$Line <- as.factor(data$Line)
data$Time <- as.factor(data$Time)
data$Image <- as.factor(data$Image)

###########################################################################
# Some simple figures for a quick overview:
###########################################################################
myPlotMeans.topleft(data$Mean, data$Line, data$Time, xlab="", ylab="", 
                    legend.lab="Time", cex=1)
myPlotMeans.topleft(data$Mean, data$Time, data$Line, xlab="", ylab="", 
                    legend.lab="Line", cex=1)
myPlotMeans.topleft(data$Area, data$Line, data$Time, xlab="", ylab="", 
                    legend.lab="Line", cex=1)

###########################################################################
# Univariate plots, sorted by mean:
###########################################################################
par(mfrow=c(2,1))
new_ord1 <- reorder(data$Line, data$Mean, FUN=mean)
bargraph.CI(new_ord1, data$Mean, ylab="Cell Size (um^2)", 
            main="06:00h Sample", las=3)

new_ord2 <- reorder(data$Line, data$Area, FUN=mean)
bargraph.CI(data$Line, data$Area, ylab="Cluster Size (um^2)", 
            main="", las=3)
