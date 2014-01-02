# ChlamyPhototaxisAnalysis.R
# Code for testing and plotting Chlamydomonas phototaxis.

# Copyright (C) 2013 K. Likos and J. Malcom

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

setwd("~/Dropbox/Chlamy_project/Chlamy wt fitness assays/wt_phototaxis_assays_Nov2012")

library(sciplot)

#FIRST, filter out all raw data with a maxiumum value
# of less than 20. 
dat <- read.table("chlamy_crosses.csv", sep=",", header=TRUE)
names(dat)
dat$ID <- as.factor(dat$ID)
dat$MEDIA <- as.factor(dat$MEDIA)

d1 <- subset.data.frame(dat, dat$MEDIA == "NEW"
                        & dat$FAR >= 20
                        | dat$MID >= 20
                        | dat$CLOSE >= 20)
#d2 <- subset.data.frame(dat, dat$MEDIA != "NEW")

#THEN, take all of the survivors and round ALL values 
# (CLOSE, MID, and FAR,) to the nearest ten. We did this  
# by making a function that does special things with signif

signif.rev <- function(x) {
  if(x < 5) {
    return(0)
  }
  if(x < 15)
    return(10)
  if(x < 100)
    return(signif(x, 1))
  if(x < 1000)
    return(signif(x, 2))
  return(signif(x, 3))
}

round_FAR <- apply(as.array(d1$FAR), MARGIN=1, FUN=signif.rev)
round_MID <- apply(as.array(d1$MID), MARGIN=1, FUN=signif.rev)
round_CLOSE <- apply(as.array(d1$CLOSE), MARGIN=1, FUN=signif.rev)
round_FAR
round_MID
round_CLOSE

PI <- (round_CLOSE - round_FAR)/(round_CLOSE + round_MID + round_FAR)

PI
ID2 <- as.character(d1$ID)
ID2
dat2 <- data.frame(ID2, PI)
dat2

dat2$PI <- as.numeric(dat2$PI)

#NOW, hopefully, we can graph this edited data...


new <- reorder(dat2$ID2, dat2$PI, FUN=mean)
bargraph.CI(new, dat2$PI, ylim=c(-1, 1), xlab="Line", 
            ylab="phototactic index", las=3)
box()
abline(h=c(1, 0.5, 0.25, 0, -0.25, -0.5), col="gray")

hist(d1$INDEX, breaks=24)

# let's check if there is a correlation between variance and intensity
names(d1)

sum_intensity <- d1$FAR + d1$MID + d1$CLOSE
d1 <- cbind(d1, sum_intensity)
plot(d1$INDEX ~ d1$sum_intensity)
PI_vars <- tapply(d1$INDEX, INDEX=d1$ID, FUN=var)
PI_mean <- tapply(d1$INDEX, INDEX=d1$ID, FUN=mean)
mean_intn <- tapply(d1$sum_intensity, INDEX=d1$ID, FUN=mean)
plot(PI_vars/PI_mean ~ mean_intn)

#TRYING AGAIN- without rounded data- only filtered


new <- reorder(d1$ID, d1$INDEX, FUN=mean)
bargraph.CI(new, d1$INDEX, ylim=c(-1, 1), xlab="Line", 
            ylab="phototactic index", las=3)
box()
abline(h=c(1, 0.5, 0.25, 0, -0.25, -0.5), col="gray")

par(mfrow=c(1,2))
hist(dat2$PI, breaks=24)
hist(d1$INDEX, breaks=24)
