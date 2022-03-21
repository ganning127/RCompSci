# Ganning Xu
# 3/18/22
# data analysis

# Setup
rm(list=ls()) # cleans out the old stuff
setwd("~/programming/programming_ncssm/rcompsci/data_analysis") # set the current working directory
source("myfunctions.R") # import but for our own files
# read in some data
data <- read.csv("protein.csv", header=TRUE)

# diagnostics and overview
dim(data) # shows dimensions (prints rows and cols)
head(data) # shows the first 6 rows
tail(data) # shows the last 6 rows
summary(data) # shows descriptive statistics about each columns (including how much missing data there is)

# check for cleanliness
clean <- ifelse(
  complete.cases(data) == TRUE, # if a row has no missing data
  1, # value if true
  0 # value if false
)


clean # prints out clean
paste("There are ", dim(data)[1]-sum(clean), " rows with missing data") # dim(data[1]) takes the total number of rows. sum(clean) takes the number of rows with complete data.

# replcae missing data with their mean/median
cols <- c(3:11) # select the 3rd to 11 columns (this is numeric data)
for (i in cols) {
  data[, i][is.na(data[, i])] <- mean(data[, i], na.rm=TRUE)  # [rows, cols]. check if the data is missing
}

# fix missing values in character (string) data (in col 2)
# we have to impute the data (infer data that is realistic)
data <- random.impute.data.frame(data, c(2)) # create new col of imputed. for multiple cols can do c(1,2,3:6)
data$Import <- data$Import.imputed # replace the original data col with the inputed one created above ($ accesses a column in table)
data$Import.imputed <- NULL # remove the last column b/c we already copied it over

# we have clean data now
library(lubridate)
data$Date <- mdy(data$Date) # mdy is a lubridate function

# do stuff to analyze data
# all data that you're using assume that you're normally distributed, so you'll need to normally distribute it if it isn't
hist(data$RedMeat) # not normally distributed
hist(log10(data$RedMeat)) # a little bit more normally distributed

hist(log10(data$WhiteMeat)) # log here doesn't really work well. the whitemean is still not normally distributed
hist(rz.transform(data$WhiteMeat)) # rz came from myfunctions.R
hist(rz.transform(data$RedMeat)) 
data[cols] <- lapply(data[cols], rz.transform) # cols was defined before. lapplly applies some math function to all functions that were named

# graph the data
library(plotluck)
plotluck(data, .~1)
