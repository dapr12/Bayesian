#----------------------------------------------
# R code for Bayesian Computation Webinar
#------------------------

#----------------------------------------------
# PART II: Posterior-Predictive Model Checking
# Require packages pscl, arm, MASS
#----------------------------------------------

require(pscl)
data(bioChemists)

# data for Poisson log linear model

fem <- as.numeric(bioChemists$fem=="Women")
mar <- as.numeric(bioChemists$mar=="Married")
kid5 <- bioChemists$kid5
phd <- bioChemists$phd
ment <- bioChemists$ment
y <- bioChemists$art

# use function bayesglm in arm package to fit a Poisson
# log-linear model with vague prior on beta

library(arm)
fit <- bayesglm(y ~ fem + mar + kid5 + phd + ment,
                family="poisson")

# simulate 1000 draws from posterior of beta

S <- sim(fit, n.sims = 1000)

# use one simulate beta from posterior to
# simulate one predictive sample ys using same covariates
# compute the standard deviation of the sample ys

post.pred.sim <- function(j){
  X <- cbind(1, fem, mar, kid5, phd, ment)
  lambda <- exp(X %*% coef(S)[j, ])
  ys <- rpois(length(y), lambda)
  sd(ys)
}

# repeat for all simulated draws of beta, store sd's

post.pred.sds <- sapply(1:1000, post.pred.sim)

# construct histogram of posterior predictive sd's
# overlay sd of observed y -- demonstrates overdispersion

library(MASS)
truehist(post.pred.sds, xlim=c(1.25, 2))
abline(v=sd(y), lwd=3, col="red")
text(1.85, 5, "OBSERVED", col="red")
