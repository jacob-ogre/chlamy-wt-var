base <- "/Users/jacobmalcom/Dropbox/Chlamy_project/Chlamy_wt_fitness_assays/extracted_data/"
d4 <- read.table(paste(base, "GxE_blups_mod4.tab", sep=""),
                  sep="\t", header=TRUE)
tail(d4)
# d4 <- d4[,-1]
# names(d4)
# d4 <- d4[,-7]
# dim(d4)

d3 <- read.table(paste(base, "GxE_blups_mod3.tab", sep=""),
                  sep="\t", header=TRUE)
tail(d3)

plot(d4$X.Intercept. ~ d3$X.Intercept.)
