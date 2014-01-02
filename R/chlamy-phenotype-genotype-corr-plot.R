# chlamy-phenotype-genotype-corr-plot.R
# Plot phenotypic in lower triangle, genotypic in upper
# Copyright (C) 2014 Kyle Hernandez

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

# NOTE:
# DO NOT CHANGE THE ORDER OF ANY OF THE INPUTS FOR MY CORRGRAM SOURCE
# THINGS BE HARDCODED YO!

base <- "~/Dropbox/Chlamy_project/Chlamy_wt_fitness_assays/extracted_data/"
source("corrgram-adjusted.r")
cor.dat <- read.table(paste(base, "Chlamy_wt_max.csv", sep=""),
                      sep=",",
                      header=TRUE,
                      check.names=FALSE)
gcor.dat <- read.csv(paste(base, "Chlamy_wt_max.csv", sep=""),
                     sep=",",
                     header=TRUE,
                     check.names=FALSE)

pretty.palette <- colorRampPalette(c("darkgreen", "white", "navy"), 
                                   space = "rgb")

# NEED TO ADD LEGEND SOMEHOW...TODO manually
corrgram(cor.dat[,3:length(cor.dat)], 
         gcor.dat[,2:length(gcor.dat)],
         lower.panel=panel.shade,
         upper.panel=panel.shade, 
         col.regions=pretty.palette, cex.labels = 0.4)
