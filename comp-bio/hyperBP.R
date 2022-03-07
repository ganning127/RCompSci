# R script for analyzing QTL data
# YOUR NAME
# DATE
# Description of script
# clean things up, set the working directory
rm(list=ls())
setwd("~/programming/programming_ncssm/rcompsci/comp-bio")
#
# load the QTL library.  NOTE!  this must be INSTALLED first in the command line, using
#  install.packages("qtl").  Once it's installed, you do not have to install again.  You do,
# however, need to tell the script that you need it for this script to work!
library(qtl)

#
# Load the data!
cross <- read.cross("csv", file="hyper.csv", genotypes=c("A", "H", "B"), na.strings="-", alleles=c("A", "B"))
jittermap(cross)
# take a look at my data, make sure it's pretty clean
cross <- est.rf(cross)
plotRF(cross) # plot the glowing line thing
plotMap(cross) # plot genetic map
plotMissing(cross) # shows that there are a lot of missing data in hyper.csv
bp <- cross$pheno$bp
hist(bp) # if this isn't normally distributed, using hist(log(bp)) can help
qqnorm(bp)
qqline(bp, col='red') # there's a good line of best fit, so it means that the data is pretty normally distributed. If it isn't you go fix it

#
#  Set up for the main scan by creating probability datasets 
# and simulated datasets (both genotype probabilities)
#  Documentation:  https://www.rdocumentation.org/packages/qtl/versions/1.50/topics/calc.genoprob
#  Documentation:  https://www.rdocumentation.org/packages/qtl/versions/1.50/topics/sim.geno
cross <- calc.genoprob(cross, step=2.0, off.end=0, error.prob=1.0e-4, map.function="haldane",stepwidth = 'fixed')
cross <- sim.geno(cross, step=2.0, off.end=0, error.prob=1.0e-4, map.function="haldane",stepwidth = 'fixed')


#
# Perform the mainscan for the QTL
cross.scanBP <- scanone(cross, pheno.col=1, model="normal", method="em")
cross.scanBP.perm <- scanone(cross, pheno.col=1, model="normal", method="em", n.perm=100) # the higher the better (typically 1000 or more)

#
# plot the mainscan
plot(cross.scanBP, main="Mainscan plot of BP") # shows that chromosone 1 is involved
thresh <- summary(cross.scanBP.perm, alpha=c(0.37, 0.1, 0.005)) #63, 90, 95% confidence intervals
abline(h=thresh[1],col='blue')
abline(h=thresh[2],col='purple')
abline(h=thresh[3],col='green')


#
#  Take a look at the summary of what you have found
summary(cross.scanBP, perm=cross.scanBP.perm, alpha=0.05) # units are in centimorgans

#
# do an effect plot.  We have two markers, so we'll do two
first <- find.marker(cross, chr=1, pos=49.2)
effectplot(cross, pheno.col=1, mname1=first) # homozogous dominant mice have higher blood pressure

second <- find.marker(cross, chr=4, pos=29.5)
effectplot(cross, pheno.col=1, mname1=second)

#
#  Create a confidence interval map, looking to zoom in on where genes might be
CIChr1 <- bayesint(cross.scanBP, lodcolumn=1, chr=1, prob=0.95) # 95 confidence interval
plot(cross.scanBP, chr=1, main='Confidence Interval for Chr1')
lines(x=CIChr1[c(1,3), 2], y=c(0, 0), type='l', col='green', lwd=4) # if you want to look for the gene, look for in green area
eCIChr1 <- bayesint(cross.scanBP, lodcolumn=2, chr=1, prob=0.95) # 95 confidence interval


CIChr4 <- bayesint(cross.scanBP, lodcolumn=1, chr=4, prob=0.95) # 95 confidence interval
plot(cross.scanBP, chr=4, main='Confidence Interval for Chr4')
lines(x=CIChr4[c(1,3), 2], y=c(0, 0), type='l', col='green', lwd=4) # if you want to look for the gene, look for in green area
eCIChr4 <- bayesint(cross.scanBP, lodcolumn=2, chr=4, prob=0.95) # 95 confidence interval
CIChr4[c(1,3), 2]

#
#  BIC MODELING
#  Is there any causality?  Does one factor, like sex, affect the phenotype?  Use BIC
# models (Baysian Information Criterion), which uses a linear regression (y=mx+b, lm) to
# see if there is causality:  Cause --> Effect
# Form:   
#  BIC(lm(effect~1))  establishes a baseline for the effect.  
#  BIC(lm(effect~cause))
#  if the second number is 10 or more LESS than the baseline number, then there is a
#  statistical evidence of causality. 



# EOF (end of file)


