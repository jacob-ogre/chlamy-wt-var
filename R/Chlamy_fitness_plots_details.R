# PlotFitnessData.R
# Plots fitness data from 96-well plate experiments.
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


library("lubridate")
library("sciplot")
library("corrgram")
source("myPlotMean.general.R")
source("myPlotMean.general.nolegend.R")
setwd(paste("~/Dropbox/Chlamy_project/Chlamy_wt_fitness_assays/", 
            "wt_fitness_assays_Oct2012/extracted_data/", sep=""))

###############################################################################
## load data
###############################################################################
dat <- read.table("Chlamy_wt_fit_exp_Oct2012.csv", header=T, sep=",")
dat <- dat[order(dat[,1], dat[,4], dat[,3]),]
dat$date <- as.factor(dmy(dat$date))
dat$plate <- as.factor(dat$plate)
dat$line <- as.factor(dat$line)
dat$treatment <- as.factor(dat$treatment)
dat$read_type <- as.factor(dat$read_type)

names(dat)
levels(dat$line)

###############################################################################
## Plot all lines in all environments (mean population all days)             ##
###############################################################################
outbas <- "~/Dropbox/writings/mss/Chlamydomonas/Wild-type_variation/Figures/"
outfil <- paste(outbas, "Figure5.pdf", sep="")

pdf(outfil, height=15, width=18)

par(mfrow=c(5,6), mar=c(3,3,1.5,1))
for(i in levels(dat$treatment)) {
    if(i != "N++100") {
        sub <- subset.data.frame(dat, dat$treatment == i)
        cur.line <- reorder(sub$line, sub$measure, FUN=mean)
        bargraph.CI(cur.line, sub$measure, xlab="", ylab="",
                    cex.names=0.7, las=2, main=i)
    }
}
rm(i, sub, cur.line)

dev.off()

###############################################################################
## How big are fitness differences (not OD)?                                 ##
###############################################################################
TAPmin <- subset.data.frame(dat, dat$treatment == "TAPmin")
tapply(TAPmin$density, TAPmin$line, FUN=mean)
rm(TAPmin)
# BIG

###############################################################################
## Extract maximum population size for each line in a replicate, plot        ##
###############################################################################
reps <- 24 * length(levels(dat$treatment)) * 2
treatment <- rep(NA, reps)       # treatment of current rep
repl <- rep(NA, reps)       # plate (proxy of rep number) of current rep
line <- rep(NA, reps)       # line of current rep
max.measure <- rep(NA, reps)       # max for current rep

count = 1
for (t in levels(dat$treatment)) {
  sub1 <- subset.data.frame(dat, dat$treatment==t)
  for (p in levels(dat$plate)) {
    sub2 <- subset.data.frame(sub1, sub1$plate==p)
    for (l in levels(dat$line)) {
      sub3 <- subset.data.frame(sub2, sub2$line==l)
      if (length(sub3$date) > 0) {
        cur.max <- max(sub3$measure)
        treatment[count] <- t
        repl[count] <- p
        line[count] <- l
        max.measure[count] <- cur.max
        count <- count + 1
      }
    }
  }
}
maxes <- data.frame(treatment, repl, line, max.measure)
write.csv(maxes, file="Chlamy_wt_assay_max_ODs.csv")

par(mfrow=c(5,7), mar=c(3,4,1.5,1))
for (t in levels(dat$treatment)) {
  sub <- subset.data.frame(maxes, maxes$treatment == t)
  cur.line <- reorder(sub$line, sub$max.measure, FUN=mean)
  bargraph.CI(cur.line, sub$max.measure, ylim=c(0, 1200), 
              xlab="", ylab="", cex.names=0.7, las=2, main=t)
}

rm(sub1, sub2, sub3, cur.max, count,reps, t, p, l, line, max.measure,
   repl, treatment)

###############################################################################
## Should also look at histograms                                            ##
###############################################################################

par(mfrow=c(4,6))
for (l in levels(dat$line)) {
  sub <- subset.data.frame(dat, dat$line == l)
  hist(sub$measure, main=l)
}

par(mfrow=c(5,7))
for (t in levels(dat$treatment)) {
  sub <- subset.data.frame(dat, dat$treatment == t)
  hist(sub$measure, main=t)
}
rm(sub, t, l, cur.line)
par(mfrow=c(1,1))

###############################################################################
## How about some phenotypic and genetic correlations?                       ##
###############################################################################
cor.dat <- read.table("Chlamy_wt_max.csv", sep=",", header=TRUE)
gcor.dat <- read.table("Chlamy_wt_max_linemeans.csv", sep=",", header=TRUE)
# remove N++100, which is synonymous with N++
gcor.dat <- gcor.dat[,-7]

# corrgram(cor.dat[,3:length(cor.dat)], order=TRUE, lower.panel=panel.shade,
#          upper.panel=panel.pie)
# corrgram(cor.dat[,3:length(cor.dat)], order=TRUE, lower.panel=panel.ellipse,
#          upper.panel=panel.pie)
# 
# corrgram(gcor.dat[,2:length(gcor.dat)], order=TRUE, lower.panel=panel.shade,
#          upper.panel=panel.pie)
# corrgram(gcor.dat[,2:length(gcor.dat)], order=TRUE, lower.panel=panel.shade,
#          upper.panel=panel.pie, abs=TRUE)

###########################################################################
# Some misc plots...
###########################################################################
par(mfrow=c(2,2), mar=c(4.5,4.5,2,1))
plot(cor.dat$P50, cor.dat$N50, xlab="P50", ylab="N50")
plot(cor.dat$Zn300, cor.dat$S150, xlab="Zn300", ylab="S150")
plot(cor.dat$NaCl_8, cor.dat$N..70, xlab="NaCl 8mg/L", ylab="N++50%")
plot(cor.dat$Zn300, cor.dat$N10, xlab="Zn300", ylab="N10")
# plot(cor.dat$S150, cor.dat$N10, xlab="S150", ylab="N10")

###############################################################################
## Gradients...                                                              ##
###############################################################################
max.dat <- read.table("Chlamy_wt_max.csv", sep=",", header=TRUE)
max.2 <- read.table("Chlamy_wt_assay_max_ODs.csv", sep=",", header=TRUE)
names(max.dat)
names(max.2)

max.2$treatment <- as.factor(max.2$treatment)
max.2$line <- as.factor(max.2$line)

####---N gradient!---####
N.grad.dat <- data.frame()
N10 <- subset.data.frame(max.2, max.2$treatment == "N10")
N50 <- subset.data.frame(max.2, max.2$treatment == "N50")
N100 <- subset.data.frame(max.2, max.2$treatment == "N++")
N150 <- subset.data.frame(max.2, max.2$treatment == "N150")

N.grad.dat <- rbind(N10, N50, N100, N150)
N.grad.dat <- droplevels(N.grad.dat)
N.grad.dat$treatment = factor(N.grad.dat$treatment,
                              levels(N.grad.dat$treatment)[c(2,4,1,3)])
rm(N10, N50, N100, N150)

myPlotMean.general(N.grad.dat$max.measure, N.grad.dat$treatment, 
                   N.grad.dat$line, xlab="Treatment", ylab="OD 440/680 (max)",
                   legend.lab="line", cex=1)

N.mapp.dat <- subset.data.frame(N.grad.dat, N.grad.dat$line=="1690" | 
  N.grad.dat$line=="2343" | N.grad.dat$line=="2938" | N.grad.dat$line=="2290")
N.mapp.dat <- droplevels(N.mapp.dat)

myPlotMean.general(N.mapp.dat$max.measure, N.mapp.dat$treatment, 
                   N.mapp.dat$line, xlab="Treatment", ylab="OD 440/680 (max)",
                   legend.lab="line", cex=1, pos="bottomright")

####---P gradient!---####
P10 <- subset.data.frame(max.2, max.2$treatment == "P10")
P50 <- subset.data.frame(max.2, max.2$treatment == "P50")
P100 <- subset.data.frame(max.2, max.2$treatment == "N++")
P150 <- subset.data.frame(max.2, max.2$treatment == "P150")
P.grad.dat <- rbind(P10, P50, P100, P150)
P.grad.dat <- droplevels(P.grad.dat)
levels(P.grad.dat$treatment)
P.grad.dat$treatment = factor(P.grad.dat$treatment,
                              levels(P.grad.dat$treatment)[c(2,4,1,3)])
rm(P10, P50, P100, P150)

P.mapp.dat <- subset.data.frame(P.grad.dat, P.grad.dat$line=="1690" | 
  P.grad.dat$line=="2343" | P.grad.dat$line=="2938" | P.grad.dat$line=="2290")
P.mapp.dat <- droplevels(P.mapp.dat)

myPlotMean.general(P.mapp.dat$max.measure, P.mapp.dat$treatment, 
                   P.mapp.dat$line, xlab="Treatment", ylab="OD 440/680 (max)",
                   legend.lab="line", cex=1, pos="topleft")
myPlotMean.general.nolegend(P.grad.dat$max.measure, P.grad.dat$treatment, 
                   P.grad.dat$line, xlab="Treatment", ylab="OD 440/680 (max)",
                   legend.lab="line", cex=1, pos="topleft")

####---S gradient!---####
S10 <- subset.data.frame(max.2, max.2$treatment == "S10")
S50 <- subset.data.frame(max.2, max.2$treatment == "S50")
S100 <- subset.data.frame(max.2, max.2$treatment == "N++")
S150 <- subset.data.frame(max.2, max.2$treatment == "S150")
S.grad.dat <- rbind(S10, S50, S100, S150)
S.grad.dat <- droplevels(S.grad.dat)
levels(S.grad.dat$treatment)
S.grad.dat$treatment = factor(S.grad.dat$treatment,
                              levels(S.grad.dat$treatment)[c(2,4,1,3)])
rm(S10, S50, S100, S150)

S.mapp.dat <- subset.data.frame(S.grad.dat, S.grad.dat$line=="1690" | 
  S.grad.dat$line=="2343" | S.grad.dat$line=="2938" | S.grad.dat$line=="2290")
S.mapp.dat <- droplevels(S.mapp.dat)

myPlotMean.general(S.mapp.dat$max.measure, S.mapp.dat$treatment, 
                   S.mapp.dat$line, xlab="Treatment", ylab="OD 440/680 (max)",
                   legend.lab="line", cex=1, pos="topleft")

####---S gradient alternative (mean pop)!---####
# S10 <- subset.data.frame(dat, dat$treatment == "S10")
# S50 <- subset.data.frame(dat, dat$treatment == "S50")
# S100 <- subset.data.frame(dat, dat$treatment == "N++")
# S150 <- subset.data.frame(dat, dat$treatment == "S150")
# S.grad.dat <- rbind(S10, S50, S100, S150)
# S.grad.dat <- droplevels(S.grad.dat)
# levels(S.grad.dat$treatment)
# S.grad.dat$treatment = factor(S.grad.dat$treatment,
#                               levels(S.grad.dat$treatment)[c(2,4,1,3)])
# rm(S10, S50, S100, S150)
# 
# S.mapp.dat <- subset.data.frame(S.grad.dat, S.grad.dat$line=="1690" | 
#   S.grad.dat$line=="2343" | S.grad.dat$line=="2938" | S.grad.dat$line=="2290")
# S.mapp.dat <- droplevels(S.mapp.dat)
# 
# myPlotMean.general(S.mapp.dat$measure, S.mapp.dat$treatment, 
#                    S.mapp.dat$line, xlab="Treatment", ylab="OD 440/680 (max)",
#                    legend.lab="line", cex=1, pos="topleft")

####---pH gradient!---####
levels(max.2$treatment)
pH55 <- subset.data.frame(max.2, max.2$treatment == "pH5.5")
pH73 <- subset.data.frame(max.2, max.2$treatment == "pH7.3")
pH91 <- subset.data.frame(max.2, max.2$treatment == "pH9.1")
pH.grad.dat <- rbind(pH55, pH73, pH91)
pH.grad.dat <- droplevels(pH.grad.dat)
levels(pH.grad.dat$treatment)
rm(pH55, pH73, pH91)

pH.mapp.dat <- subset.data.frame(pH.grad.dat, pH.grad.dat$line=="1690" | 
  pH.grad.dat$line=="2343" | pH.grad.dat$line=="2938" | 
  pH.grad.dat$line=="2290")
pH.mapp.dat <- droplevels(pH.mapp.dat)

myPlotMean.general(pH.mapp.dat$max.measure, pH.mapp.dat$treatment, 
                   pH.mapp.dat$line, xlab="Treatment", ylab="OD 440/680 (max)",
                   legend.lab="line", cex=1, pos="bottomleft", ylim=c(450, 1000))
myPlotMean.general.nolegend(pH.grad.dat$max.measure, pH.grad.dat$treatment, 
                   pH.grad.dat$line, xlab="Treatment", ylab="OD 440/680 (max)",
                   legend.lab="line", cex=1, pos="bottomleft", ylim=c(450, 1000))

par(mfrow=c(2,2))
myPlotMean.general.nolegend(N.mapp.dat$max.measure, N.mapp.dat$treatment, 
                   N.mapp.dat$line, xlab="", ylab="OD 440/680 (max)",
                   legend.lab="line", cex=1, pos="bottomright")
myPlotMean.general.nolegend(P.mapp.dat$max.measure, P.mapp.dat$treatment, 
                   P.mapp.dat$line, xlab="", ylab="",
                   legend.lab="line", cex=1, pos="topleft")
myPlotMean.general.nolegend(S.mapp.dat$max.measure, S.mapp.dat$treatment, 
                   S.mapp.dat$line, xlab="Treatment", ylab="OD 440/680 (max)",
                   legend.lab="line", cex=1, pos="topleft")
myPlotMean.general(pH.mapp.dat$max.measure, pH.mapp.dat$treatment, 
                   pH.mapp.dat$line, xlab="Treatment", ylab="",
                   legend.lab="line", cex=1, pos="bottomleft")

