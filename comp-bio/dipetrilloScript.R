# R script for analyzing QTL data
# Ganning Xu
# 3/6/22
# clean things up, set the working directory
rm(list=ls())
setwd("~/programming/programming_ncssm/rcompsci/comp-bio")

library(qtl) # import the qtl library


cross <- read.cross("csv", file="dipetrillo.csv", genotypes=c("C", "H", "S"), na.strings="-", alleles=c("C", "S"))
jittermap(cross)
cross <- est.rf(cross)
plotRF(cross) # plot the glowing line thing
plotMap(cross) # plot genetic map
plotMissing(cross) # shows that there are a lot of missing data in hyper.csv


# BLOOD PRESSURE
bp <- cross$pheno$Mean.BP
hist(bp) # if this isn't normally distributed, using hist(log(bp)) can help
qqnorm(bp)
qqline(bp, col='red') 
cross <- calc.genoprob(cross, step=2.0, off.end=0, error.prob=1.0e-4, map.function="haldane",stepwidth = 'fixed')
cross <- sim.geno(cross, step=2.0, off.end=0, error.prob=1.0e-4, map.function="haldane",stepwidth = 'fixed')
cross.scanBP <- scanone(cross, pheno.col=3, model="normal", method="em")
cross.scanBP.perm <- scanone(cross, pheno.col=3, model="normal", method="em", n.perm=100) # the higher the better (typically 1000 or more)
plot(cross.scanBP, main="Mainscan plot of BP") 
thresh <- summary(cross.scanBP.perm, alpha=c(0.10, 0.05, 0.01)) # 90, 95, 99% confidence intervals
abline(h=thresh[1],col='blue')
abline(h=thresh[2],col='purple')
abline(h=thresh[3],col='green')
summary(cross.scanBP, perm=cross.scanBP.perm, alpha=0.05) # units are in centimorgans

first <- find.marker(cross, chr=1, pos=66.5) 
effectplot(cross, pheno.col=3, mname1=first) # create an 1st effect plot
second <- find.marker(cross, chr=16, pos=48.3)
effectplot(cross, pheno.col=3, mname1=second) # create 2nd effect plot

CIChr1 <- bayesint(cross.scanBP, lodcolumn=1, chr=1, prob=0.95) # 95 confidence interval
plot(cross.scanBP, chr=1, main='Confidence Interval for Chr1')
lines(x=CIChr1[c(1,3), 2], y=c(0, 0), type='l', col='green', lwd=4) # if you want to look for the gene, look for in green area
eCIChr1 <- bayesint(cross.scanBP, lodcolumn=1, chr=1, prob=0.95) # 95 confidence interval

CIChr16 <- bayesint(cross.scanBP, lodcolumn=1, chr=16, prob=0.95) # 95 confidence interval
plot(cross.scanBP, chr=16, main='Confidence Interval for Chr16')
lines(x=CIChr16[c(1,3), 2], y=c(0, 0), type='l', col='green', lwd=4) # if you want to look for the gene, look for in green area
CIChr16 <- bayesint(cross.scanBP, lodcolumn=1, chr=16, prob=0.95) # 95 confidence interval
CIChr16[c(1,3), 2]

# HEART RATE
hr <- cross$pheno$Mean.HR
hist(hr) # if this isn't normally distributed, using hist(log(bp)) can help
cross <- calc.genoprob(cross, step=2.0, off.end=0, error.prob=1.0e-4, map.function="haldane",stepwidth = 'fixed')
cross <- sim.geno(cross, step=2.0, off.end=0, error.prob=1.0e-4, map.function="haldane",stepwidth = 'fixed')
cross.scanHR <- scanone(cross, pheno.col=4, model="normal", method="em")
cross.scanHR.perm <- scanone(cross, pheno.col=4, model="normal", method="em", n.perm=100) # the higher the better (typically 1000 or more)
plot(cross.scanHR, main="Mainscan plot of HR") 
thresh <- summary(cross.scanHR.perm, alpha=c(0.10, 0.05, 0.01)) # 90, 95, 99% confidence intervals
abline(h=thresh[1],col='blue')
abline(h=thresh[2],col='purple')
abline(h=thresh[3],col='green')
summary(cross.scanHR, perm=cross.scanHR.perm, alpha=0.05) # units are in centimorgans

# HEART WEIGHT
wt <- cross$pheno$Heart.Wt
hist(wt) # if this isn't normally distributed, using hist(log(bp)) can help
