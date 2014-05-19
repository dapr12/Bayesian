# Data: wine ratings, wine prices, and review words from http://www.tastings.com

wine = read.table("wine.tbl", header = T)

# Columns 3 and up are predictors (words occurring in reviews).
# So we have (num columns) - 2 predictors

predictor.col = 3:ncol(wine)
num.predictors = ncol(wine) - 2

# add a column PriceClass that cuts up prices into three bins:
# cheap, medium, and expensive

wine$PriceClass = "medium"
wine[wine$WinePrice < 10,]$PriceClass = "cheap"
wine[wine$WinePrice > 25,]$PriceClass = "expensive"
wine$PriceClass = as.factor(wine$PriceClass)

# Do Naive Bayes classification.

library(klaR)
nb = NaiveBayes(wine[,predictor.col], wine$PriceClass)

# exploring the nb object:
# "tables" has the probabilities
# P(predictor = TRUE | class = "cheap") 
# (and likewise for the other target classes).
# "apriori" has the probabilities
# P(class = "cheap") and likewise for the other target classes.

names(nb)
nb$apriori
nb$tables
nb$tables$toasty
length(nb$tables)

# make a data frame that has one row for each predictor
# with columns "predictor",
# "cheap" with P(predictor = TRUE | class = "cheap"),
# "expensive" with P(predictor = TRUE | class = "expensive")

cheap.prob.v = c()
for(i in 1:num.predictors) cheap.prob.v = append(cheap.prob.v, nb$tables[[i]][1,2])

expensive.prob.v = c()
for(i in 1:num.predictors) expensive.prob.v = append(expensive.prob.v, nb$tables[[i]][2,2])

word.prob = data.frame(predictor = names(nb$tables), cheap = cheap.prob.v, expensive = expensive.prob.v)

# what are the words w (predictors) with the highest P(w = T |cheap) and P(w = T|expensive)?

head(word.prob[order(word.prob$cheap, decreasing=T),], n=20)
head(word.prob[order(word.prob$expensive, decreasing=T),], n=20)

# looking at the ratio of cheap versus expensive probability.

word.prob$cheap.exp.ratio = word.prob$cheap / word.prob$expensive

# these are the words where the expensive probability is highest compared to the cheap probability

head(word.prob[order(word.prob$cheap.exp.ratio),], n=20)

# these are the words where the cheap probability is highest compared to the expensive probability

head(word.prob[order(word.prob$cheap.exp.ratio,decreasing=T),], n=20)

# separating training and test data
# for some reason I could not get this to work with categorical predictors
# so I convert T/F to 1/0

wine2 = wine
for (i in predictor.col) wine2[,i] = as.numeric(wine2[,i])

# Use the last 9/10 of the data for training, and the first 1/10 for test
nrow(wine)
nrow(wine)/10

nb2 = NaiveBayes(wine2[168:nrow(wine2),predictor.col], wine2[168:nrow(wine2),]$PriceClass)

# prediction for the first wine

predict(nb2, wine2[1, predictor.col])

# what did it actually cost?

wine2[1,c("WinePrice", "PriceClass")]

# predictions for all the test data

predictions = predict(nb2, wine2[1:167,predictor.col])
names(predictions)

# posterior probabilities: always with one very strong strongest class 
# Probabilities for the "winning" class are almost always 
# between 85% and 99%. 
# This is typical for Naive Bayes.

predictions$posterior

# view the prediction and the true price next to each other

cbind(as.character(predictions$class), as.character(wine2[1:167,]$PriceClass))

# and look at the original texts

wine.desc = read.table("wine-orig.tbl", header=T)
head(wine.desc)

# Can you guess the price from the descriptions?
# In particular look at wine 24, which was mis-classified as expensive when in fact it was cheap
