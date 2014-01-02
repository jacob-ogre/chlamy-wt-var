# Chlamy_GxE_updated.R
# Analysis of Chlamy wt GxE.
# Written in Aug 2013 by Kyle Hernandez, kmhernan@utexas.edu
# Updated Jan 2014 by Jacob Malcom, jmalcom@uconn.edu
# Copyright (C) 2013 K. Hernandez

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


library(lme4)

###########################################################################
# load data
###########################################################################
base <- "~/Dropbox/Chlamy_project/Chlamy_wt_fitness_assays/"
datf <- paste(base, "extracted_data/wt_fitness_maxF_FINAL.tab", sep="")
dat <- read.table(datf, header=TRUE, sep="\t")

summary(dat)
names(dat)
keep <- c("1373", "1690", "2246", "2290", "2342", "2343", "2344", "2607",
          "2608", "2931", "2932", "2935", "2936", "2937", "2938", "4414",
          "89", "90")
dat <- droplevels(subset(dat, dat$Line %in% keep))
dat <- droplevels(subset(dat, dat$Env != "NULL"))
dat$Coord <- as.factor(dat$Coord)
dat$Env <- as.factor(dat$Env)

###########################################################################
# Analysis
###########################################################################
# I think this is right -- KMH
# I also think this is right -- JWM
fit.1 <- lmer(dat$Fluor ~ 1 + (1|dat$Line))
fit.2 <- lmer(dat$Fluor ~ 1 + dat$Env + (1|dat$Line))
fit.3 <- lmer(dat$Fluor ~ 1 + dat$Env + (1|dat$Line) + (1|dat$Line:dat$Env))

# add plate coords for update:
fit_4 <- lmer(dat$Fluor ~ 1 + dat$Env + (1|dat$Coord) + (1|dat$Line) + 
              (1|dat$Line:dat$Env))

fit_aov <- anova(fit.1, fit.2, fit.3, fit_4)
fit_aov
summary(fit.3)
summary(fit_4)

ran.3 <- ranef(fit.3)
tab.gxe <- ran.3$`dat$Line:dat$Env`

fix_4 <- fixef(fit_4)
ran_4 <- ranef(fit_4)
new_gxe <- ran_4$`dat$Line:dat$Env`

plot(new_gxe$`(Intercept)`, tab.gxe$`(Intercept)`)
cor.test(new_gxe$`(Intercept)`, tab.gxe$`(Intercept)`)

# Check correlations
# library(reshape)
# max.dat <- read.table("~/Dropbox/Chlamy_project/Chlamy_wt_fitness_assays/extracted_data/Chlamy_wt_max_linemeans.csv", check.names=FALSE, header = T, sep = ",")
# blup.dat <- read.table("~/Dropbox/Chlamy_project/Chlamy_wt_fitness_assays/extracted_data/formatted_gxe_blups.tab", header = F)
# names(blup.dat) <- c("line", "treatment", "value")

# new.max.dat <- melt(max.dat, id=c("line"))
# names(new.max.dat) <- c("line", "treatment", "value")

# combined.dat <- merge(new.max.dat, blup.dat, by=c("line", "treatment"))
# cor(combined.dat$value.x, combined.dat$value.y)
# par(mfrow=c(1,2))
# plot(abs(combined.dat$value.y) ~ combined.dat$value.x,  ylab="abs(Blups)", xlab="Line Means", main="Absolute")
# plot(combined.dat$value.y ~ combined.dat$value.x,  ylab="Blups", xlab="Line Means", main="Raw")

# Analysis of individual lines for stats for Fig 6
N_list <- c("N++", "N10", "N50", "N150")
N_grad <- droplevels(subset(dat, dat$Env %in% N_list))
fit.1 <- lmer(N_grad$Fluor ~ 1 + (1|N_grad$Line))
fit_N <- lmer(N_grad$Fluor ~ 1 + (1|N_grad$Coord) + (1|N_grad$Line) +
              (1|N_grad$Line:N_grad$Env))
fit_aov <- anova(fit.1, fit_N)
fit_aov
ranef(fit_N)
c2936 <- subset(N_grad, N_grad$Line == "2936")
mod2936 <- lm(c2936$Fluor ~ c2936$Env)
summary(mod2936)

P_list <- c("N++", "P10", "P50", "P150")
P_grad <- droplevels(subset(dat, dat$Env %in% P_list))
fit.1 <- lmer(P_grad$Fluor ~ 1 + (1|P_grad$Line))
fit_P <- lmer(P_grad$Fluor ~ 1 + (1|P_grad$Coord) + (1|P_grad$Line) +
              (1|P_grad$Line:P_grad$Env))
fit_aov <- anova(fit.1, fit_P)
fit_aov
ranef(fit_P)
c2290 <- subset(P_grad, P_grad$Line == "2290")
mod2290 <- lm(c2290$Fluor ~ c2290$Env)
summary(mod2290)

S_list <- c("N++", "S10", "S50", "S150")
S_grad <- droplevels(subset(dat, dat$Env %in% S_list))
fit.1 <- lmer(S_grad$Fluor ~ 1 + (1|S_grad$Line))
fit_S <- lmer(S_grad$Fluor ~ 1 + (1|S_grad$Coord) + (1|S_grad$Line) +
              (1|S_grad$Line:S_grad$Env))
fit_aov <- anova(fit.1, fit_S)
fit_aov
ranef(fit_S)
c2932 <- subset(S_grad, S_grad$Line == "2932")
mod2932 <- lm(c2932$Fluor ~ c2932$Env)
summary(mod2932)

c2937 <- subset(S_grad, S_grad$Line == "2937")
mod2937 <- lm(c2937$Fluor ~ c2937$Env)
summary(mod2937)

###########################################################################
# Write results to file
###########################################################################
outf <- paste(base, "extracted_data/GxE_blups_rev.tab", sep="")
write.table(tab.gxe, outf, sep = "\t", quote=FALSE, row.names=TRUE)
