# myPlotMean.general.R
# A generic function to plot a response variable by two factors.
#
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

# This is my function to create a plot of means without Rcmdr, on which 
# this code is based.
myPlotMean.general <- function(response, fact1, fact2, xlab, ylab, legend.lab, 
                               cex, ylim=c(0,1000), pch=1:(1+nlevs2), 
                               lty=1:nlevs2, pos="topleft", error.width=NULL, 
                               par.mar=NULL) {
	xlab
	ylab
	legend.lab
	cex
    pos
    if (is.null(error.width)) {
        error.width <- 0.125
    }
    if (is.null(par.mar)) {
        par.mar=c(5, 8, 3, 2)
    }
	means <- tapply(response,list(fact1,fact2),mean)
	sds <- tapply(response,list(fact1,fact2),sd)
	ns <- length(response)
	sds <- sds/sqrt(ns)
	yrange <- c(min(means-sds,na.rm=TRUE), max(means+sds,na.rm=TRUE))
	levs1 <- levels(fact1)
	levs2 <- levels(fact2)
	nlevs1 <- length(levs1)
	nlevs2 <- length(levs2)
	if (length(lty) == 1) {
      	lty <- rep(lty, nlevs2)
    }
	par(cex.lab=cex, cex.axis=cex, font.lab=2, mar=par.mar, lwd=1)
	plot(c(1, nlevs1),yrange,type="n",xlab="", ylab="",axes=FALSE)
	title(xlab=xlab)
	title(ylab=ylab, line=3)
	box()
	axis(2)
	axis(1, at=1:nlevs1, labels=levs1, las=3)
	for (i in 1:nlevs2) {
		points(1:nlevs1, means[,i], type="b", pch=pch[i], cex=cex, lty=lty[i])
		arrows(1:nlevs1, means[,i]+(1.96*sds[,i]), 1:nlevs1, 
               means[,i]-(1.96*sds[,i]), angle=90, code=3, lty=lty[i], 
               length=error.width) 
    }
	x.posn <- 1
	y.posn <- 1
	legend(pos, levs2, pch=pch, lty=lty, title=legend.lab, cex=cex-0.1)
}
