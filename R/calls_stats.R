# Plot count data. */
# Copyright (C) 2013 Kyle Hernandez */
#
# This program is free software; you can redistribute it and/or modify */
# it under the terms of the GNU General Public License as published by */
# the Free Software Foundation; either version 2 of the License, or */
# (at your option) any later version. */
#
# This program is distributed in the hope that it will be useful, */
# but WITHOUT ANY WARRANTY; without even the implied warranty of */
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the */
# GNU General Public License for more details. */
#
# You should have received a copy of the GNU General Public License */
# along with this program; if not, see <http://www.gnu.org/licenses/>. */


library(lattice)

data <- read.delim("/Users/kmhernan/Desktop/chlamy/new/data_count.txt", 
                   header = TRUE, sep = "\t")
barchart(ratio ~ sample | enzyme, data = data,
	scales = list(x=list(cex=0.75, rot=45)),
	ylab = "Het/Total Ratio")

barchart(homs ~ sample | enzyme, data = data,
	scales = list(x=list(cex=0.75, rot=45)),
	ylab = "Hom Counts")
	
barchart(hets ~ sample | enzyme, data = data,
	scales = list(x=list(cex=0.75, rot=45)),
	ylab = "Het Counts")
