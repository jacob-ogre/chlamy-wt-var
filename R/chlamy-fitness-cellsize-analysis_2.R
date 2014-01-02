# Testing for phenotype-fitness relationships.
# Copyright (C) 2013 K. Hernandez and J. Malcom

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

base <- "~/Dropbox/Chlamy_project/"
cellSZ <- read.table(paste(base, "/Chlamy_photos/Chlamy_2_time_point_all.csv",
                           sep=""), 
                     sep=",", header = TRUE)

base2 <- paste(base, 
               "Chlamy_wt_fitness_assays/wt_fitness_assays_Oct2012/extracted_data/", 
               sep="")
fitness <- read.table(paste(base2, "Chlamy_wt_assay_max_ODs.csv", sep=""),
                      sep=",", header=TRUE)

cellSZ$Line <- as.factor(cellSZ$Line)
cellSZ$Time <- as.factor(cellSZ$Time)
cellSZ$Image <- as.factor(cellSZ$Image)

npp <- fitness[ which(fitness$treatment=="N++" &
                      fitness$line !="33A5" &
                      fitness$line !="33C1" &
                      fitness$line !="33D4" &
                      fitness$line !="43B2" &
                      fitness$line !="43B6" &
                      fitness$line !="43C5"),]
npp <- droplevels(npp)
npp$line <- factor(npp$line, 
                   levels(npp$line)[c(17, 18, 1, 2, 3, 4, 5, 6, 7, 8, 
                                      9, 10, 11, 12, 13, 14, 15, 16)],
                   order=T)
npp$treatment <- as.factor(npp$treatment)

dat <- aggregate(npp$max.measure, by=list(npp$line), mean)
dat <- cbind(dat, aggregate(cellSZ$Mean, by=list(cellSZ$Line), mean)[2])
names(dat) <- c("Line", "Fitness", "Size")

plot(dat$Fitness ~ dat$Size)
fm.1 <- lm(dat$Fitness ~ dat$Size)
anova(fm.1)

hist(dat$Size, 30)

