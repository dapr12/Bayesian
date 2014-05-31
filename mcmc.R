#
# An example of MCMC inference with R.
#
# Time-stamp: <2010-06-12 18:04 Petri Koistinen>

# Run this code a few lines at a time using, e.g., the R editor.

#------------------
# The data
#------------------

y <- c(-0.91, 1.27, 0.85, 0.89, -0.22, 1.07, 2.25, -1.20, -0.95, -0.82)
n <- length(y)

#-----------------------------------------------------------------
# Our model states that
#
# [Y[i] | mu, psi] i.i.d~ Normal(mu, 1 / psi),   i = 1, ..., n
#
# Parameter vector theta = (mu, psi), where mu is the expected value
# and psi is the precision of the normal population.
#
# We use the improper prior with density proportional to 1/psi.
#-----------------------------------------------------------------


# Calculate crude guesses for the parameters

print(mu.crude <- mean(y))
print(sigma.crude <- sd(y))
print(psi.crude <- 1 / var(y))

#--------------------------------------------------------------------
# Define functions for calculating  log-likelihood, log-prior, and
# log(unnormalized posterior)
#--------------------------------------------------------------------

# Log-likelihood
loglik <- function(theta) {
  sum(dnorm(log = TRUE, y, mean = theta[1], sd = 1 / sqrt(theta[2])))
}

# log-prior for the improper prior
logprior <- function(theta) (-log(theta[2]))

# Log unnormalized posterior
logpost <- function(theta) (loglik(theta) + logprior(theta))

# By trial and error, I have set up a grid in the parameter space,
# which covers most of the posterior probability mass.

mu.grid <- seq(-2.5, 2.5, len = 101)
psi.grid <- seq(0.05, 2.6, len = 201)

# Now we calculate the value of the unnormalized posterior
# on the grid.

poste <- matrix(0, nrow = length(mu.grid), ncol = length(psi.grid))
for (i in seq(along = mu.grid))
  for (j in seq(along = psi.grid))
    poste[i, j] <- exp(logpost(c(mu.grid[i], psi.grid[j])))

# This gives some idea of the shape of the posterior

contour(mu.grid, psi.grid, poste)

# However, I want to specify better levels at which to draw
# the contour lines.

maxposte <- max(poste)
# Contours are drawn at the following levels; notice that we define
# the levels relative to the maximum of the (positive!) function
levs <- c(0.9, 0.5, 0.1, 0.01) * maxposte

# This makes a nice plot of the posterior.
contour(mu.grid, psi.grid, poste, levels = levs, drawlabels = FALSE)

# ... Hmm, interesting. Markedly non-normal shape.

#-----------------------------------
# The Gibbs sampler
#-----------------------------------

# Number of iterations
niter <- 2000

# Allocate a (niter + 1) x 2 matrix for the simulations.

sim.gibbs <- matrix(0, nrow = niter + 1, ncol = 2)

# Set meaningful names to the columns of the matrix sim.gibbs:
colnames(sim.gibbs) <- c('mu', 'psi')

# Initial values
mu.cur <- mu.crude + rnorm(1, sd = 0.2)
psi.cur <- psi.crude * exp(rnorm(1, sd = 0.2))
sim.gibbs[1, 1] <- mu.cur
sim.gibbs[1, 2] <- psi.cur

ybar <- mean(y)
for (i in 0:niter) {
  mu.cur <- rnorm(1, mean = ybar, sd = 1 / sqrt(n * psi.cur))
  psi.cur <- rgamma(1, n/2, 0.5 * sum( (y - mu.cur)^2 ))
  sim.gibbs[i + 1, 1] <- mu.cur
  sim.gibbs[i + 1, 2] <- psi.cur
}

#----------------------
# Visualize results
#----------------------

# Make the contour plot and add the simulated points.
contour(mu.grid, psi.grid, poste, levels = levs, drawlabels = FALSE)
points(sim.gibbs[, 'mu'], sim.gibbs[, 'psi'], pch = '.')
title(main = 'Gibbs sampler', xlab = 'mu', ylab = 'psi')


#--------------------------------------
# Calculate posterior summaries:
#--------------------------------------

# Actually, I am not really interested in (mu, psi) but the parameters
# (mu, sigma), where sigma = 1 / sqrt(psi)

# We get samples for these parameters as follows:

mu.sig.gibbs <- cbind(sim.gibbs[, 'mu'], 1 / sqrt(sim.gibbs[, 'psi']))

# ... cbind() makes a matrix out of column vectors

# Set meaningful names to the columns of matrix mu.sig.gibbs:

colnames(mu.sig.gibbs) <- c('mu', 'sigma')

# Let's reject first half of the simulations

mu.sig.gibbs <- mu.sig.gibbs[-(1 : round(niter/2) ), ]

# Posterior means:

apply(mu.sig.gibbs, 2, mean)

# This is a fancy way of calculating
# mean(mu.sig.gibbs[,1]) and mean(mu.sig.gibbs[,2])
# with a single call.

# Posterior covariance:

cov(mu.sig.gibbs)

# Quantiles of the marginal posterior distributions:

apply(mu.sig.gibbs, 2, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))

# You can obtain the same information with two calls:
# quantile(mu.sig.gibbs[,1], probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
# quantile(mu.sig.gibbs[,2], probs = c(0.025, 0.25, 0.5, 0.75, 0.975))


#---------------------------------------------------------
#
# Reparametrize the model in terms of 
# phi = (mu, lambda), where lambda = log(psi)
#
#---------------------------------------------------------

# Define functions for calculating the log-likelihood,
# and log(unnormalized posterior) after the reparametrization.


r.loglik <- function(phi) loglik(c(phi[1], exp(phi[2])))
r.logpost <- function(phi) 
  (logpost(c(phi[1], exp(phi[2]))) + phi[2])


# Setup a grid for the lambda-axis
lambda.grid = seq(-2.5, 1.5, len = 101)

# Evaluate the unnormalized log-posterior on the grid.

r.poste <- matrix(0, nrow = length(mu.grid), ncol = length(lambda.grid))
for (i in seq(along = mu.grid))
  for (j in seq(along = lambda.grid))
    r.poste[i, j] <- exp(r.logpost(c(mu.grid[i], lambda.grid[j])))

maxposte <- max(r.poste)
r.levs <- c(0.9, 0.5, 0.1, 0.01) * maxposte

# This makes a nice plot of the posterior in the reparametrized model:


contour(mu.grid, lambda.grid, r.poste, levels = r.levs, 
drawlabels = FALSE)

# Clearly non-normal even after the reparametrization.
# However, let's try to do a normal approximation:

#---------------------------------------
# Normal approximation to the posterior
#---------------------------------------

# Search the posterior mode in the reparametrized model.
# Ask the optimizer optim() to return the Hessian matrix
# which is found by numerical differentiation.

?optim

# We use the crude guesses as initial values.
# We search for the minimizer of the negative log-posterior.
# Neldear-Mead is a method which does not need a function
# for calculating the gradient of the object function.

neg.r.logpost <- function(phi) (-r.logpost(phi))


r.postmode <- optim(c(mu.crude, log(psi.crude)), neg.r.logpost, 
method = "Nelder-Mead", hessian = TRUE)

# This is the posterior mode (in the reparametrized model)

r.postmode$par

# This is the covariance matrix of the normal approximation
# to the posterior (in the reparametrized model)
# Function solve calculates the inverse matrix.

r.V <- solve(r.postmode$hessian)
rownames(r.V) <- c('mu', 'lambda')
colnames(r.V) <- c('mu', 'lambda')
print(r.V)


#------------------------------------------------------------
# Metropolis-Hastings sampler for the reparametrized model
#------------------------------------------------------------



sim.mh <- matrix(0, nrow = niter + 1, ncol = 2)
colnames(sim.mh) <- c('mu', 'lambda')
mu.cur <- mu.crude + rnorm(1, sd = 0.5)
lambda.cur <- log(psi.crude) + rnorm(1, sd = 0.5)

# I use a two-variate normal random walk proposal N(0, S), where the
# covariance of the proposal distribution is proportional
# to the covariance matrix of the normal approximation
# to the posterior.

A <- t(chol(r.V))

# Now A %*% t(A) is equal to the covariance matrix C
# of the normal approximation to the posterior.
# I can draw a sample from N(0, b^2 * C) as follows,
#
# b * A %*% rnorm(2) 
#
# Matrix multiplication is denoted by %*% and t() means the transpose.
# R interprets vectors as row or column vectors -- 
# depending on the context.

# Initial values.

phi.cur <- c(mu.crude, log(psi.crude)) + rnorm(2, sd = 0.1)
sim.mh[1, ] = phi.cur

# Keep track of the number of accepted proposals
nacc <- 0
b <- 2

for (i in 0:niter) {
  phi.prop <- phi.cur + b * A %*% rnorm(2)

  log.r <- r.logpost(phi.prop) - r.logpost(phi.cur)

  if (runif(1) < exp(log.r)) {
    nacc <- nacc + 1
    phi.cur <- phi.prop
  }
  sim.mh[i + 1, ] <- phi.cur
}

# ... my initial guess for b, b = 2, gives a reasonable
# acceptance rate.  No need for further tuning.

print(nacc / niter)

#--------------------
# Visualize results
#--------------------

contour(mu.grid, lambda.grid, r.poste, levels = r.levs,
drawlabels = FALSE)
points(sim.mh[, 'mu'], sim.mh[, 'lambda'], pch = '.')
title('MH sampler', xlab = 'mu', ylab = 'log(psi)')

#-------------------------------
# Calculate posterior summaries
#-------------------------------

# Check whether the posterior mean is about the same as the
# numerically identified posterior mode, and whether the
# posterior covariance of (mu, lambda) is about the same
# as the covariance matrix of the normal approximation.

print(r.postmode$par)
apply(sim.mh, 2, mean)

print(r.V)
cov(sim.mh)

# ... the posterior is so non-normal that we should not
# expect very good agreement.

# Actually, I am interested in the parameters
# (mu, sigma), where sigma = exp(-lambda / 2)

mu.sig.mh <- cbind(mu = sim.mh[, 'mu'], sigma = exp( - sim.mh[, 'lambda'] / 2))

# Let's reject first half of the iterations

mu.sig.mh <- mu.sig.mh[-(1 : round(niter/2) ), ]

# Posterior means

apply(mu.sig.mh, 2, mean)

# Posterior covariance

cov(mu.sig.mh)

# Quantiles of the marginal posterior distributions

apply(mu.sig.mh, 2, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))




#------------------------------------------------
# MCMC output analysis with the coda library
# which you have to install first ...
#------------------------------------------------

library(coda)

# This library has an awkward menu interface

# codamenu()

# convert the arrays to objects of type mcmc

gibbs.mcmc <- mcmc(mu.sig.gibbs)
summary(gibbs.mcmc)
plot(gibbs.mcmc)
autocorr.plot(gibbs.mcmc)
effectiveSize(gibbs.mcmc)

mh.mcmc <- mcmc(mu.sig.mh)
summary(mh.mcmc)
plot(mh.mcmc)
autocorr.plot(mh.mcmc)
effectiveSize(mh.mcmc)

# My original Gibbs sampler is clearly better than my M--H sampler.

#------------------------------------------------
# Bayesian output analysis with the boa library,
# which you have to install first ...
#------------------------------------------------

library(boa)

# The library has an awkward menu interface.

# boa.menu()

# Luckily, we can access its functions also directly:

help(package = "boa")

boa.init()
boa.chain.add(boa.importMatrix('mu.sig.gibbs'), 'Gibbs')
boa.chain.add(boa.importMatrix('mu.sig.mh'), 'MH1')
boa.chain.add(boa.importMatrix('mu.sig.mh2'), 'MH2')
boa.plot('trace')
boa.plot('acf')

# Hmm ... the Gibbs sampler seems to be the best.
# I would choose it for serious analysis of the data.

# Summary statistics about the chains:

boa.print.stats()

# The column 'MC Error' shows an estimate of the Monte Carlo
# standard error for the 'Mean' column.  

boa.quit()
 
