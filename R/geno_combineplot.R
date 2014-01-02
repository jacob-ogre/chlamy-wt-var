# Plot multiple marker sets.
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


# read in datasets from multiple markersets (i.e., RAD and SNP markers)
# lengths.tab is a text file of chromosome lengths
dg<-read.delim("../datasets/par_conv_new_genotype.tab", header=T)
dl<-read.delim("../datasets/lengths.tab", header=F)
odg<-read.delim("../datasets/cat_oldmarker.csv", header = T)
chrl<-levels(factor(dg$Chr))

# Sort on Chromosome and Position
dg<- dg[order(dg$Chr, dg$Position),]
odg<- odg[order(odg$Chr, odg$Position),]
pdf("../figures/brachy_combine_pretty.pdf", width = 8.5, height = 11)
par(mar=c(0.05, 2.15, 0.05, 0.01))

# set your rows based on how many samples you have
# It will plot each dataset in a single row
par(mfrow=c(152, 1))
par(cex=0.5)

# The chromosomes are named different things in my dataset here
# because i didn't make the SNP marker dataset
# Currently, It only does one chromosome/scaffold at a time, but you could change it

dc <- subset(dg,Chr=="Chr1")
odc<- subset(odg,Chr=="1")

# for samples in one dataset
for (a in 9:159){
	
	# for samples in the other dataset
	for (b in 4:168){
		# Match sample IDs (each dataset has to have the same naming
		# scheme for individuals)
		if(colnames(dc[a])==colnames(odc[b])){
			plot(0,0,xlim=c(0,dl[1,2]),
				ylim=c(-1,1), ylab="", xlab="", type="n",
				fg=0, col.axis=0)
				mtext(colnames(odc)[b], line=0, side=2, las=1, xpd=NA, cex=0.4)
			# this will plot grey behind everything so N's will be
			# grey			
			rect(0, -.9,dl[1,2], 0.9, border = "NA", col="grey70", cex=1)
			
			# Here, I loop through each row of one dataset
                	for (c in 1:nrow(odc[b])){
                        	if (c!=nrow(odc[b])){
					# I skip the highest position cause there isn't c+1

					# Now if the current position and the next position
					# are equal plot a rect based on the type
                                	if (odc[c,b]=="a" && odc[c+1,b]=="a"){
                                        	rect(odc$Position[c], -.9, odc$Position[c+1], 0.9,
                                                	border="NA",col="red", cex=1)
					} else if (odc[c,b]=="b" && odc[c+1,b]=="b"){
                                        	rect(odc$Position[c], -.9, odc$Position[c+1], 0.9,
                                                	border="NA", col="blue", cex=1)
					} else if(odc[c,b]=="h" && odc[c+1,b]=="h"){
                                        	rect(odc$Position[c], -0.9, odc$Position[c+1], 0.9,
                                                	border="NA", col="green", cex=1)
					}
                        	} 
				# I guess this is me testing a way to fill in the final marker
				else if (c==nrow(odc[b])){
					  if (odc[c-1,b]=="a" && odc[c,b]=="a"){
                                                rect(odc$Position[c-1], -.9, odc$Position[c], 0.9,
                                                        border="NA",col="red", cex=1)
                                        } else if (odc[c-1,b]=="b" && odc[c,b]=="b"){
                                                rect(odc$Position[c-1], -.9, odc$Position[c], 0.9,
                                                        border="NA", col="blue", cex=1)
                                        } else if(odc[c-1,b]=="h" && odc[c,b]=="h"){
                                                rect(odc$Position[c-1], -0.9, odc$Position[c], 0.9,
                                                        border="NA", col="green", cex=1)
                                        } 
				}
			}
			
			# do same for next dataset
			plot(0,0,xlim=c(0,dl[1,2]),
                        	ylim=c(-1,1), ylab="", xlab="", type="n",
                               	fg=0, col.axis=0)
			mtext(colnames(dc)[a], line=0, side=2, las=1, xpd=NA, cex=0.4)
			rect(0, -.9,dl[1,2], 0.9, border = "NA", col="grey70", cex=1)
			for (d in 1:nrow(dc[a])){
                                if (d!=nrow(dc[a])){
                                       	if (dc[d,a]=="A" && dc[d+1,a]=="A"){
                                               	rect(dc$Position[d], -.9, dc$Position[d+1],
							 0.9, border="NA",col="red", cex=1)
                                       	} else if (dc[d,a]=="B" && dc[d+1,a]=="B"){
                                               	rect(dc$Position[d], -.9, dc$Position[d+1],
 							0.9, border="NA", col="blue", cex=1)
                                       	} else if(dc[d,a]=="HET" && dc[d+1,a]=="HET"){
                                               	rect(dc$Position[d], -0.9, dc$Position[d+1], 
							0.9, border="NA", col="green", cex=1)
                                       	}
                                } else if (d==nrow(odc[a])){
					  if (odc[d-1,a]=="A" && odc[d,a]=="A"){
                                                rect(odc$Position[d-1], -.9, odc$Position[d], 0.9,
                                                        border="NA",col="red", cex=1)
                                        } else if (odc[d-1,a]=="B" && odc[d,a]=="B"){
                                                rect(odc$Position[d-1], -.9, odc$Position[d], 0.9,
                                                        border="NA", col="blue", cex=1)
                                        } else if(odc[d-1,a]=="HET" && odc[c,a]=="HET"){
                                                rect(odc$Position[d-1], -0.9, odc$Position[d], 0.9,
                                                        border="NA", col="green", cex=1)
                                        } 
				}
			} 
		plot.new()
        	}        
        }
}
dev.off()
