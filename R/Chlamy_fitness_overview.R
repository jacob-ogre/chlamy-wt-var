# Chlamy_fitness_overview.R
# Plot a general overview of Chlamy fitness data.
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


source("./myPlotMedian.general.R")

###########################################################################
# load the data
###########################################################################
file <- ####!--FIND FILE PATH IF THIS IS NEEDED!########
dat <- read.table(file, sep=",", header=T)
as.Date(dat$date, "%d/%m/%y")
is.factor(dat$date)
names(dat)

fluor <- subset.data.frame(dat, dat$read_type=="fluor")
abs <- subset.data.frame(dat, dat$read_type=="abs")

a1690 <- subset(abs, abs$line=="1690")
a2938 <- subset(abs, abs$line=="2938")
a2931 <- subset(abs, abs$line=="2931")
f1690 <- subset(fluor, fluor$line=="1690")
f2938 <- subset(fluor, fluor$line=="2938")
f2931 <- subset(fluor, fluor$line=="2931")

d1690 <- 0.4718 * exp(f1690$measure * 0.007)
ad1690 <- 0.6798 * exp(a1690$measure * 2.9819)
d2938 <- 0.4718 * exp(f2938$measure * 0.007)
d2931 <- 0.4718 * exp(f2931$measure * 0.007)

###########################################################################
# Plot the data
###########################################################################
myPlotMedian.general(a1690$measure, a1690$date, a1690$treatment,
                     xlab="", ylab="OD (750nm Abs)", legend.lab="Treatment",
                     cex=1, pos="topright")

myPlotMedian.general(ad1690, a1690$date, a1690$treatment,
                     xlab="", ylab="Density (cells per 0.625ul)", 
                     legend.lab="Treatment", cex=1, pos="topright")

myPlotMedian.general(a2938$measure, a2938$date, a2938$treatment,
                     xlab="", ylab="OD (750nm Abs)", legend.lab="Treatment",
                     cex=1, pos="topright")

myPlotMedian.general(d1690, f1690$date, f1690$treatment,
                     xlab="", ylab="OD (750nm Abs)", legend.lab="Treatment",
                     cex=1, pos="topright")

myPlotMedian.general(d2938, f2938$date, f2938$treatment,
                     xlab="", ylab="OD (750nm Abs)", legend.lab="Treatment",
                     cex=1, pos="topright")

myPlotMedian.general(d2931, f2931$date, f2931$treatment,
                     xlab="", ylab="OD (750nm Abs)", legend.lab="Treatment",
                     cex=1, pos="topright")
