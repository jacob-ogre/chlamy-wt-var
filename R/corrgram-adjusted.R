# corrgram.r
# Time-stamp: <06 Nov 2012 15:01:49 c:/x/rpack/corrgram/R/corrgram.r>

# Author: Kevin Wright
# Copyright 2006-2012 Kevin Wright
# License: GPL2

# The corrgram function was derived from the 'pairs' function.
# Code for plotting ellipses was derived from the ellipse package.
# Additional ideas from the plot.corr function in the 'arm' package.

# To do: Add a legend/ribbon

##################################################
### Adjusted to work for our pheno/geno corr plot
### by K Hernandez 2013
##################################################

#
# This really only works for two non-corr matrices that are equal size.
# I have it hard coded to order based on the 'y' matrix (genetic corr in our case)
# which is also expected to be the 'top' panel.

corrgram <-
  function (x, y, type=NULL,
            labels, panel = panel.shade, #...,
            lower.panel = panel, upper.panel = panel,
            diag.panel = NULL, text.panel = textPanel,
            label.pos = 0.5, label.srt=0,
            cex.labels = NULL, font.labels = 1,
            row1attop = TRUE, dir="", gap = 0,
            abs=FALSE,
            col.regions = colorRampPalette(c("red","salmon","white","royalblue","navy")),
            ...) {
  
  if(is.null(order)) order <- FALSE
  
  # Direction
  if(dir=="") {
    if(row1attop) dir <- "left" else dir <- "right"
  }
  if (dir=="\\") dir <-  "left"
  if (dir=="/") dir <-  "right"  
  
  if (ncol(x) < 2 || ncol(y) < 2) stop("Only one column in the argument to 'corrgram'")
  
  # Do we have a data.frame or correlation matrix?
  maybeCorr <- FALSE
  
  if(is.null(type)){
    if(maybeCorr)
      type <- "corr"
    else
      type <- "data"
  } else if(type=="data"){
    if(maybeCorr)
      warning('This looks like a correlation matrix.')
  } else if(type=="cor" | type=="corr") {
    type <- "corr"
    if(!maybeCorr)
      stop('This is NOT a correlation matrix.')
  } else {
    stop("unknown data type in 'corrgram'")
  }

  # Remove non-numeric columns from data frames
  if(type=="data" & !is.matrix(x)) x <- x[ , sapply(x, is.numeric)]
  if(type=="data" & !is.matrix(y)) y <- y[ , sapply(y, is.numeric)]
  
  # If a data matrix, then calculate the correlation matrix
  xmat <- cor(x, use="pairwise.complete.obs")
  ymat <- cor(y, use="pairwise.complete.obs")
  xmat <- if(abs) abs(xmat) else xmat
  ymat <- if(abs) abs(ymat) else ymat
  
  # Order by angle size between PCAs (first two) of correlation matrix
  ### KH - ORDER ONLY BY 'Y' DATASET

  y.eigen <- eigen(ymat)$vectors[,1:2]
  e1 <- y.eigen[,1]
  e2 <- y.eigen[,2]
  alpha <- ifelse(e1>0, atan(e2/e1), atan(e2/e1)+pi)
  ord <- order(alpha)
  y <- if(type=="data") y[,ord] else y[ord, ord]


  textPanel <- function(x = 0.5, y = 0.5, txt, cex, font, srt) {
    text(x, y, txt, cex=cex, font=font, srt=srt)
  }
  
  localAxis <- function(side, x, y, xpd, bg, col=NULL, main, oma, ...) {
    ## Explicitly ignore any color argument passed in as
    ## it was most likely meant for the data points and
    ## not for the axis.
    if(side %%2 == 1) Axis(x, side=side, xpd=NA, ...)
    else Axis(y, side=side, xpd=NA, ...)
  }
  
  # Don't pass some arguments on to the panel functions via the '...'
  localPlot <- function(..., main, oma, font.main, cex.main)
    plot(...)
  localLowerPanel <- function(..., main, oma, font.main, cex.main)
    lower.panel(...)
  localUpperPanel <- function(..., main, oma, font.main, cex.main)
    upper.panel(...)
  localDiagPanel <- function(..., main, oma, font.main, cex.main)
    diag.panel(...)

  dots <- list(...)
  nmdots <- names(dots)

  # Check for non-numeric data
  if (!is.matrix(x)) {
    x <- as.data.frame(x)
    for(i in seq(along=names(x))) {
      if(is.factor(x[[i]]) || is.logical(x[[i]]))
        x[[i]] <- as.numeric(x[[i]])
      if(!is.numeric(unclass(x[[i]])))
        stop("non-numeric argument to 'corrgram'")
    }
  } else if (!is.numeric(x)) stop("non-numeric argument to 'corrgram'")

  ### KH - CHECK Y MATRIX
  if (!is.matrix(y)) {
    y <- as.data.frame(y)
    for(i in seq(along=names(y))) {
      if(is.factor(y[[i]]) || is.logical(y[[i]]))
        y[[i]] <- as.numeric(y[[i]])
      if(!is.numeric(unclass(y[[i]])))
        stop("non-numeric argument to 'corrgram'")
    }
  } else if (!is.numeric(y)) stop("non-numeric argument to 'corrgram'")
  # Get panel functions
  panel <- match.fun(panel)
  if((has.lower <- !is.null(lower.panel)) && !missing(lower.panel))
    lower.panel <- match.fun(lower.panel)
  if((has.upper <- !is.null(upper.panel)) && !missing(upper.panel))
    upper.panel <- match.fun(upper.panel)
  
  has.diag  <- !is.null(diag.panel)
  if(has.diag && !missing( diag.panel))
    diag.panel <- match.fun( diag.panel)
  
  if(dir=="left") {
    tmp <- lower.panel; lower.panel <- upper.panel; upper.panel <- tmp
    tmp <- has.lower; has.lower <- has.upper; has.upper <- tmp
  }

  # Plot layout
  ## KH - Check if the dimensions of the matrices are the same
  nc <- ncol(x)
  nt <- ncol(y)
  if(nc != nt)
    stop("unequal matrices!")
  has.labs <- TRUE
  ## KH - dealing with the different sorting of the two matrices
  if (missing(labels)) {
    labels <- colnames(y)
    labelsx <- colnames(x)
    if (is.null(labels)) labels <- paste("var", 1:nc)
  }
  else if(is.null(labels)) has.labs <- FALSE
  if(is.null(text.panel)) has.labs <- FALSE
  
  oma <- if("oma" %in% nmdots) dots$oma else NULL
  main <- if("main" %in% nmdots) dots$main else NULL
  
  if (is.null(oma)) {
    oma <- c(4, 4, 4, 4)
    if (!is.null(main)) oma[3] <- 6 # Space for the title
  }
  opar <- par(mfrow = c(nc, nc), mar = rep.int(gap/2, 4), oma = oma)
  on.exit(par(opar))

  # Main loop to draw each panel
  for (i in if(dir=="left") 1:nc else nc:1)
    for (j in 1:nc) {
      ## KH - hardcoded crap to make sure the same traits are being plotted between
      ## the pheno matrix and the geno matrix
      yi_lab <- labels[i]
      yj_lab <- labels[j]
      xi <- match(yi_lab, labelsx)
      xj <- match(yj_lab, labelsx)
      # Set up plotting area
      ## KH - use xi/xj in place of 'i' and 'j' for localPlots
      localPlot(x[, xj], x[, xi], xlab = "", ylab = "", axes = FALSE, type = "n", ...)
      if(i == j || (i < j && has.lower) || (i > j && has.upper) ) {
        
        if(i == j) {
          # Diagonal panel
          if (has.diag) {
            if(type=="data")
              localDiagPanel(as.vector(x[, xi]), NULL, ...)
            else
              localDiagPanel(NULL, x[xi,xi], ...)
          }
          
          # Diagonal text
          if (has.labs) {
            par(usr = c(0, 1, 0, 1))
            if(is.null(cex.labels)) {
              l.wid <- strwidth(labels, "user")
              cex.labels <- max(0.8, min(2, .9 / max(l.wid)))
            }
            text.panel(0.5, label.pos, labels[i],
                       cex = cex.labels, font = font.labels, srt=label.srt)
          }
        } else if(i < j) {
          # Lower panel
          ## KH - Hardcoded to be phenotype matrix
          if(type=="data")
            localLowerPanel(as.vector(x[, xj]), as.vector(x[, xi]), NULL, col.regions, ...)
          else
            localLowerPanel(NULL, NULL, x[xj,xi], col.regions, ...)
        } else {
          # Upper panel
          ## KH - hardcoded to be genotype matrix. Note that it's ok to use 'i' and 'j' here
          if(type=="data")
            localUpperPanel(as.vector(y[, j]), as.vector(y[, i]), NULL, col.regions, ...)
          else
            localUpperPanel(NULL, NULL, y[j,i], col.regions, ...)
        }
        
      } else { # No panel drawn
        par(new = FALSE)
      }
    }
  if (!is.null(main)) {
    font.main <- if("font.main" %in% nmdots) dots$font.main else par("font.main")
    cex.main <- if("cex.main" %in% nmdots) dots$cex.main else par("cex.main")
    mtext(main, 3, 3, TRUE, 0.5, cex = cex.main, font = font.main)
  }
  invisible(NULL)
}

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# Panel functions

panel.pts <- function(x, y, corr=NULL, col.regions, ...){
  
  # For correlation matrix, do nothing
  if(!is.null(corr)) return()
  
  plot.xy(xy.coords(x, y), type="p", ...)
  box(col="lightgray")
}

panel.pie <- function(x, y, corr=NULL, col.regions, ...){

  # Coordinates of box
  usr <- par()$usr
  minx <- usr[1]; maxx <- usr[2]
  miny <- usr[3]; maxy <- usr[4]
  # Multiply the radius by .97 so the circles do not overlap
  rx <- (maxx-minx)/2 * .97
  ry <- (maxy-miny)/2 * .97
  centerx <- (minx+maxx)/2
  centery <- (miny+maxy)/2

  segments <- 60
  angles <- seq(0,2*pi,length=segments)
  circ <- cbind(centerx + cos(angles)*rx, centery + sin(angles)*ry)
  lines(circ[,1], circ[,2], col='gray30',...)

  # If corr is NULL, then we calculate it
  if(is.null(corr))
    corr <- cor(x, y, use='pair')
  
  # Overlay a colored polygon
  ncol <- 14
  pal <- col.regions(ncol)
  col.ind <- as.numeric(cut(corr, breaks=seq(from=-1, to=1, length=ncol+1),
                            include.lowest=TRUE))
  col.pie <- pal[col.ind]
  
  segments <- round(60*abs(corr),0) 
  if(segments>0){ # Watch out for the case with 0 segments
    angles <- seq(pi/2, pi/2+(2*pi* -corr), length=segments)
    circ <- cbind(centerx + cos(angles)*rx, centery + sin(angles)*ry)
    circ <- rbind(circ, c(centerx, centery), circ[1,])
    polygon(circ[,1], circ[,2], col=col.pie)
  }
  
}

panel.shade <- function(x, y, corr=NULL, col.regions, ...){
  
  if(is.null(corr))
    corr <- cor(x, y, use='pair')
  
  ncol <- 14
  pal <- col.regions(ncol)
  col.ind <- as.numeric(cut(corr, breaks=seq(from=-1, to=1, length=ncol+1),
                            include.lowest=TRUE))
  usr <- par("usr")
  # Solid fill
  rect(usr[1], usr[3], usr[2], usr[4], col=pal[col.ind], border=NA)
  # Add diagonal lines

  if(!is.na(corr)) {
    rect(usr[1], usr[3], usr[2], usr[4], density=5,
         angle=ifelse(corr>0, 45, 135), col="white")
  }
  # Boounding box needs to plot on top of the shading, so do it last.
  box(col='grey90')
}

panel.ellipse <- function(x,y, corr=NULL, col.regions, ...){

  # For correlation matrix, do nothing
  if(!is.null(corr)) return()
  
  # Draw an ellipse
  dfn <- 2
  dfd <- length(x)-1
  shape <- var(cbind(x,y),na.rm=TRUE)
  keep <- (!is.na(x) & !is.na(y))
  center <- c(mean(x[keep]),mean(y[keep]))
  radius <- sqrt(dfn*qf(.68,dfn,dfd))
  segments <- 75
  angles <- seq(0,2*pi,length=segments)
  unit.circle <- cbind(cos(angles),sin(angles))
  ellipse.pts <- t(center+radius*t(unit.circle%*%chol(shape)))
  ellx <- ellipse.pts[,1]
  elly <- ellipse.pts[,2]
  # Truncate ellipse at min/max or at bounding box
  usr <- par()$usr
  minx <- usr[1] #min(x, na.rm=TRUE)
  maxx <- usr[2] #max(x, na.rm=TRUE)
  miny <- usr[3] #min(y, na.rm=TRUE)
  maxy <- usr[4] #max(y, na.rm=TRUE)
  ellx <- ifelse(ellx < minx, minx, ellx)
  ellx <- ifelse(ellx > maxx, maxx, ellx)
  elly <- ifelse(elly < miny, miny, elly)
  elly <- ifelse(elly > maxy, maxy, elly)
  lines(ellx, elly, col='gray30',...)

  # Filled ellipse
  # polygon(ellx, elly, col="blue", ...)

  # Add a lowess line through the ellipse.  Use 'ok' to remove NAs
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) 
    lines(stats::lowess(x[ok], y[ok], f = 2/3, iter = 3), 
          col = "red", ...)  
}

panel.bar <- function(x, y, corr=NULL, col.regions, ...){
  # Use 'bars' as in Friendly, figure 1
  
  usr <- par()$usr
  minx <- usr[1]; maxx <- usr[2]
  miny <- usr[3];  maxy <- usr[4]

  if (is.null(corr)) 
    corr <- cor(x, y, use = "pair")
  ncol <- 14
  pal <- col.regions(ncol)
  col.ind <- as.numeric(cut(corr, breaks = seq(from = -1, to = 1, 
                                    length = ncol + 1), include.lowest = TRUE))
  col.bar <- pal[col.ind]
  if(corr < 0) {
    # Draw up from bottom
    maxy <- miny + (maxy-miny) *  abs(corr)
    rect(minx, miny, maxx, maxy, col = pal[col.ind], 
         border = "lightgray")
  } else if (corr > 0){
    # Draw down from top
    miny <- maxy - (maxy-miny)*corr
    rect(minx, miny, maxx, maxy, col = pal[col.ind], 
         border = "lightgray")
  } 

}

panel.conf <- function(x, y, corr=NULL, col.regions, digits=2, cex.cor, ...){

  auto <- missing(cex.cor)  
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  
  # For correlation matrix, only show the correlation
  if(!is.null(corr)) {
    est <- corr
    est <- formatC(est, digits=digits, format='f')
    if(auto) cex.cor <- 0.7/strwidth(est)
    text(0.5, 0.6, est, cex=cex.cor)
    
  } else { # Calculate correlation and confidence interval
    results <- cor.test(x, y, alternative = "two.sided")
  
    est <- results$estimate
    est <- formatC(est, digits=digits, format='f')
    if(auto) cex.cor <- 0.7/strwidth(est)
    text(0.5, 0.6, est, cex=cex.cor)
    
    ci <- results$conf.int
    ci <- formatC(ci, digits=2, format='f')
    ci <- paste("(",ci[1],",",ci[2],")",sep="")
    if(auto) cex.cor <- 0.8/strwidth(ci)
    text(0.5, 0.3, ci, cex=cex.cor)
  }
}

panel.txt <- function(x=0.5, y=0.5, txt, cex, font, srt){
  text(x, y, txt, cex=cex, font=font, srt=srt)
}

panel.density <- function(x, corr=NULL, ...){
  # For correlation matrix, do nothing  
  if(!is.null(corr)) return()

  dd = density(x, na.rm=TRUE)
  xr=range(dd$x)
  yr=range(dd$y)
  par(usr = c(min(xr), max(xr), min(yr), max(yr)*1.1))
  plot.xy(xy.coords(dd$x, dd$y), type="l", col="black", ...)
  box(col="lightgray")
}

panel.minmax <- function(x, corr=NULL, ...){
  # For correlation matrix, do nothing  
  if(!is.null(corr)) return()
  # Put the minimum in the lower-left corner and the
  # maximum in the upper-right corner
  minx <- round(min(x, na.rm=TRUE),2)
  maxx <- round(max(x, na.rm=TRUE),2)
  text(minx, minx, minx, cex=1, adj=c(0,0))
  text(maxx, maxx, maxx, cex=1, adj=c(1,1))
}

