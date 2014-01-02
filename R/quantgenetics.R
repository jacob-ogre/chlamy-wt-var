# Returns Vg, Ve, Vp, and heritability.
# Copyright (C) 2014 Jacob Malcom
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>.


quantgene <- function(x) {
    obs_mod <- lmer(INDEX ~ 1 + (1|ID), data=x)
    obs_sum <- summary(obs_mod)
    obs_slt <- slot(obs_sum, "REmat")
    obs_Vg <- as.numeric(obs_slt[5])
    obs_Ve <- as.numeric(obs_slt[6])
    obs_Vp <- obs_Vg + obs_Ve
    obs_H2 <- obs_Vg / obs_Vp
    return(list(obs_Vg, obs_Ve, obs_Vp, obs_H2))
}
