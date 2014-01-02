# Plot of GFF SNP categories.
# Copyright (C) 2013 Kyle Hernandez, kmhernan@utexas.edu

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
setwd("/home/kmhernan/Projects/Chlamy/WT_Natural_Variation/results/gff_location_tables/")

# 
# Read in data
dat.gff <- read.delim("combined_snp_location.tab", header = TRUE)
# Read in pdf outfile name
out.image <- "../figures/gff_location/gff_snp_features.pdf"

# Since we can't have fine control of axes with facet_wrap and facet_grid,
# we have to use viewports tricks.
vplayout <- function(x, y) viewport(layout.pos.row=x, layout.pos.col=y)

# Create viewports layers
# First everything but intergenic/Intronic
p1 <- ggplot(subset(dat.gff, Feature !="Intergenic" & Feature !="Intron"), aes(x=Dataset, y=Ratio, fill=Dataset)) + 
      geom_bar(stat="identity") +
      expand_limits(y=c(0,0.4)) +
      scale_fill_grey() + 
      facet_wrap(~Feature) +
      theme_bw() +
      theme(axis.text.x=element_blank(), 
            axis.ticks.x=element_blank(), 
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            legend.position="null",
            plot.margin=unit(c(0.2,0.2,0.2,0.2), "cm"))

# Now just Intergenic and Intronic
p2 <- ggplot(subset(dat.gff, Feature =="Intergenic" | Feature =="Intron"), aes(x=Dataset, y=Ratio, fill=Dataset)) + 
  geom_bar(stat="identity") +
  ylim(0, 0.75) +
  scale_fill_grey() + 
  facet_wrap(~Feature) +
  theme_bw() +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position=c(.75,.92),
        legend.key.size=unit(0.30, "cm"),
        legend.text=element_text(size=6),
        legend.title=element_blank(),
        plot.margin=unit(c(0.2,0.2,0.2,0.2), "cm"))
pdf(out.image, width=7, height=4)
grid.newpage()
pushViewport(viewport(layout=grid.layout(4,7)))
print(p1,vp=vplayout(1:4,1:4))
print(p2,vp=vplayout(1:4,5:7))
dev.off()
