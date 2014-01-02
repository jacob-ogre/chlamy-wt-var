# Plots of linkage disequilibrium.
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

library(ggplot2)
library(grid)

# Read in data
dat.LD <-  read.delim("/home/kmhernan/Projects/Chlamy/WT_Natural_Variation/results/LD_rad/CR_all_MD10000_WC250.tab", header = TRUE)
out.image <- "/home/kmhernan/Projects/Chlamy/WT_Natural_Variation/results/figures/ld_rad/ld_inset.pdf"

# Make subset dataframes
sub.1000 <- subset(dat.LD, DIST < 1000)
sub.100 <- subset(dat.LD, DIST < 100)

# Use trusty viewports function
vplayout <- function(x, y) viewport(layout.pos.row=x, layout.pos.col=y)

p1 <- ggplot(sub.1000, aes(x=DIST, y=RSQ)) +
      geom_point(size=1.5, alpha=1/2) +
      ylab(expression(paste(rho^2))) +
      xlab("Distance (bp)") +
      theme_bw() +
      theme(plot.margin=unit(c(0.2,0.2,0.2,0.2), "cm"),
            axis.text.x=element_text(size=7),
            axis.text.y=element_text(size=7),
            panel.grid=element_blank())

p2 <- ggplot(sub.100, aes(x=DIST, y=RSQ)) +
      geom_point(size=1) +
      theme_bw() + 
      theme(axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.x=element_text(size=6),
            axis.text.y=element_text(size=6),
            panel.grid=element_blank())

pdf(out.image, height=6, width=6)
grid.newpage()
pushViewport(viewport(layout=grid.layout(6,6)))
print(p1,vp=vplayout(1:6,1:6))
print(p2,vp=vplayout(2:3,3:5))
dev.off()
