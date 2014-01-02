# Plot genome-wide diversity statistics.
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

args<-commandArgs(TRUE)

# Load data
dat <- read.delim(args[1], header = TRUE)
#dat <- read.delim("/home/kmhernan/Projects/Chlamy/WT_Natural_Variation/results/pop_stats/reduced_100kb_5window.tab",
#                  header = TRUE)

ord <- c("chromosome_1", "chromosome_2", "chromosome_3", "chromosome_4",
         "chromosome_5", "chromosome_6", "chromosome_7", "chromosome_8",
         "chromosome_9", "chromosome_10", "chromosome_11", "chromosome_12",
         "chromosome_13", "chromosome_14", "chromosome_15", "chromosome_16",
         "chromosome_17")

cl <- c("black", "grey40","black", "grey40","black", "grey40","black", "grey40",
        "black", "grey40","black", "grey40","black", "grey40","black", "grey40", 
        "black")

cbg <- c("grey90", "white","grey90", "white","grey90", "white","grey90", "white",
        "grey90", "white","grey90", "white","grey90", "white","grey90", "white", 
        "grey90")

# Rounding function
roundUp <- function(x,to)
{
  to*(x%/%to + as.logical(x%%to))
}

# Set axis limits
theta_lim <- c(0, roundUp(max(dat$theta) * 1000, 10))
pi_lim <- c(0, roundUp(max(dat$pi) * 1000, 2))
s_lim <- c(0, roundUp(max(dat$S) * 100, 2))
n_chk <- min(dat$n_loci) - 500
n_lim <- c(ifelse(n_chk < 0, 0, n_chk), roundUp(max(dat$n_loci), 1000)) 

pdf(args[2], 8, 8)
par(mfrow=c(4,18), mar=c(4.5,0,0.1,0.1), oma=c(0,0,0,0))
for (i in 1:length(ord)){
	curr <- subset(dat, Chromosome == ord[i])
	if (i == 1){
	  plot.new()
	  plot(0,0, ylim=theta_lim, 
	       xlim=c(0, max(curr$Bin)+10), 
	       type="n", axes=FALSE,
	       ylab="", xlab="")
	  mtext(expression(paste("Blockwise ", theta, "w (", 10^-3, ")")), side=2, cex=.5, padj=-2.75)
	  rect(par("usr")[1], 0, par("usr")[2], theta_lim[2], col = cbg[i], border=NA)
	  lines(curr$theta * 1000 ~ curr$Bin, lwd=1.5, col=cl[i])
	  mtext(i, side=1, cex=.75)
    axis(2, cex.axis=0.75)
	} else {
	  plot(0,0, ylim=theta_lim, 
	       xlim=c(0, max(curr$Bin)+10), 
	       type="n", axes=FALSE,
	       ylab="", xlab="")
	  rect(par("usr")[1], 0, par("usr")[2], theta_lim[2], col = cbg[i], border=NA)
	  lines(curr$theta * 1000 ~ curr$Bin, lwd=1.5, col=cl[i])
	  mtext(i, side=1, cex=.75)
	}
}
# Plot pi
for (i in 1:length(ord)){
	curr <- subset(dat, Chromosome == ord[i])
	if (i == 1){
    plot.new()
	  plot(0,0, ylim=pi_lim, 
	       xlim=c(0, max(curr$Bin)), 
	       type="n", axes=FALSE,
	       ylab="", xlab="")
	  rect(par("usr")[1], 0, par("usr")[2], pi_lim[2], col = cbg[i], border=NA)
	  lines(curr$pi * 1000 ~ curr$Bin, lwd=1.5, col=cl[i])
	  mtext(i, side=1, cex=.75)
		axis(2, cex.axis=0.75)
		mtext(expression(paste("Blockwise ",pi, " (", 10^-3, ")")), side=2, cex=.5, padj=-2.75)
	} else {
	    plot(0,0, ylim=pi_lim, 
	         xlim=c(0, max(curr$Bin)), 
		      type="n", axes=FALSE,
		      ylab="", xlab="")
	    rect(par("usr")[1], 0, par("usr")[2], pi_lim[2], col = cbg[i], border=NA)
	    lines(curr$pi * 1000 ~ curr$Bin, lwd=1.5, col=cl[i])
	    mtext(i, side=1, cex=.75)
	}
}

for (i in 1:length(ord)){
  curr <- subset(dat, Chromosome == ord[i])
  if (i == 1){
    plot.new()
    plot(0,0, ylim=s_lim, 
         xlim=c(0, max(curr$Bin)), 
         type="n", axes=FALSE,
         ylab="", xlab="")
    rect(par("usr")[1], 0, par("usr")[2], s_lim[2], col = cbg[i], border=NA)
    lines(curr$S * 100 ~ curr$Bin, lwd=1.5, lty= 1, col=cl[i])
    mtext(i, side=1, cex=.75)
    axis(2, cex.axis=0.75)
    mtext("Avg. % Segregating", side=2, cex=.5, padj=-5)
  } else {
    plot(0,0, ylim=s_lim, 
         xlim=c(0, max(curr$Bin)), 
         type="n", axes=FALSE,
         ylab="", xlab="")
    rect(par("usr")[1], 0, par("usr")[2],s_lim[2], col = cbg[i], border=NA)
    lines(curr$S * 100 ~ curr$Bin, lwd=1.5, lty= 1, col=cl[i])
    mtext(i, side=1, cex=.75)
  }
}

for (i in 1:length(ord)){
  curr <- subset(dat, Chromosome == ord[i])
  if (i == 1){
    plot.new()
    plot(0,0, ylim=n_lim, 
         xlim=c(0, max(curr$Bin)), 
         type="n", axes=FALSE,
         ylab="", xlab="")
    rect(par("usr")[1], n_lim[1], par("usr")[2], n_lim[2], col = cbg[i], border=NA)
    lines(curr$n_loci ~ curr$Bin, lwd=1.5, lty=1, col=cl[i])
    mtext(i, side=1, cex=.75)
    axis(2, cex.axis=0.751)
    mtext("Avg. N Loci", side=2, cex=.5, padj=-5)
  } else {
    plot(0,0, ylim=n_lim, 
         xlim=c(0, max(curr$Bin)), 
         type="n", axes=FALSE,
         ylab="", xlab="")
    rect(par("usr")[1], n_lim[1], par("usr")[2], n_lim[2], col = cbg[i], border=NA)
    lines(curr$n_loci ~ curr$Bin, lwd=1.5, lty=1, col=cl[i])
    mtext(i, side=1, cex=.75)
  }
}
dev.off()
