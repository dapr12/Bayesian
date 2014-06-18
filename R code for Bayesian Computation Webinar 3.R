#----------------------------------------------
# R code for Bayesian Computation Webinar
# Jim Albert - June 12, 2014
#----------------------------------------------
#---------------------------------------------------------------------
# PART III: INTRODUCTION TO JAGS
# Negative binomial regression
# Example from Jackman, 2009
# Should have previously installed JAGS on the computer
# Require R packages pscl, rjags, coda, LearnBayes
#---------------------------------------------------------------------

# model description

modelString <- "
model{
for(i in 1:n){
mu[i] <- beta[1] + beta[2]*fem[i] + beta[3]*mar[i]
+ beta[4]*kid5[i] + beta[5]*phd[i]
+ beta[6]*ment[i]
lambda[i] <- exp(mu[i])
p[i] <- r/(r+lambda[i])
y[i] ~ dnegbin(p[i],r)
}

beta[1:6] ~ dmnorm(b0[1:6],B0[,])
r ~ dunif(0,50)
}"
 
# write model description to a file

writeLines(modelString, con="negbin1.bug")

# data setup

require(pscl)
data(bioChemists)

forJags <- list(n=dim(bioChemists)[1], ## sample size
                fem=as.numeric(bioChemists$fem=="Women"), ## covariates
                mar=as.numeric(bioChemists$mar=="Married"),
                kid5=bioChemists$kid5,
                phd=bioChemists$phd,
                ment=bioChemists$ment,
                y=bioChemists$art, ## response
                b0 = rep(0,6), ## prior hyperparameters
                B0 = diag(.0001,6))

# define some initial parameter values for the MCMC

inits <- list(list(beta=rep(0,6), r=1))

# run JAGS through the rjags package

require(rjags)
foo <- jags.model(file="negbin1.bug",
                  data=forJags,
                  inits=inits)

update(foo, 5000)

out <- coda.samples(foo,
                    variable.names=c("beta","r"),
                    n.iter=10000)

# summaries of all parameters

summary(out)

# look at trace plots of all parameters

xyplot(out)

##################################################
# try a different MCMC algorithm in LearnBayes
##################################################

# write a function defining the log posterior

nbreg <- function(theta, y, X){
  beta <- theta[-7]
  r <- theta[7]
  mu <- X %*% beta
  lambda <- exp(mu)
  p <- r / (r + lambda)
  sum(dnbinom(y, prob=p, size=r, log=TRUE))
}

# data setup

fem <- as.numeric(bioChemists$fem=="Women")
mar <- as.numeric(bioChemists$mar=="Married")
kid5 <- bioChemists$kid5
phd <- bioChemists$phd
ment <- bioChemists$ment
X <- cbind(1, fem, mar, kid5, phd, ment)
y <- bioChemists$art

# starting value

theta <- c(rep(0, 6), 1)

library(LearnBayes)

# Here I do three runs of laplace since this was slow to converge
# to the posterior mode starting at this initial starting value

theta <- c(rep(0, 6), 1)
fit1 <- laplace(nbreg, theta, y, X)
fit1 <- laplace(nbreg, fit1$mode, y, X)
fit1 <- laplace(nbreg, fit1$mode, y, X)

# 10,000 draws from random walk Metropolis algorithm

fit2 <- rwmetrop(nbreg, list(var=fit1$var, scale=1),
                 fit1$mode, 10000, y, X)

# we see better mixing here than JAGS

sdata <- data.frame(fit2$par)
names(sdata) <- c("beta1", "beta2", "beta3", "beta4",
                  "beta5", "beta6", "r")
summary(mcmc(sdata))
xyplot(mcmc(sdata))
