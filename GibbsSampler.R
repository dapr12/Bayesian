###################################
# Gibbs sampler for simple normal mean hierarchical model
###################################


n<-1000
mu<-3
sigma<-3
sigma2<-sigma^2
y<-rnorm(n,mu,sigma)

####################
# Prior Hyperparms #
####################
m0<-0
t0<-0     # Noninformative
a<-b<-0   # Noninformative

#########
# Inits #
#########
tau<-1
mtau<-a+n/2

#################
# Store Results #
#################
nsim<-5000
Tau<-Mu<-Sigma<-rep(0,nsim)

#########
# Gibbs #
#########
for (i in 1:nsim) {
  v<-1/(tau*n + t0)
  m<-v*(tau*sum(y)+t0*m0)
  mu<- Mu[i]<-rnorm(1,m,sqrt(v))

  tau<-rgamma(1,mtau,b+sum((y-mu)^2)/2)
  Sigma[i]<-1/tau
}
#######################
# Posterior summaries #
#######################
mean(Mu[1501:nsim])	# Discard first 1500 as burn-in
sd(Mu[1501:nsim])
mean(Sigma[1501:nsim])
sd(Sigma[1501:nsim])
par(ask=F)
plot(501:nsim,Mu[501:nsim],type="l",col="lightgreen")
abline(h=mean(Mu[501:nsim]))

##################################
# Draw mu directly from marginal #
#    posterior assuming non-     #
#      informative priors        #
##################################

mu2<-rt(3500,n-1)*sd(y)/sqrt(n)+mean(y)
mean(mu2)
sd(mu2)

par(ask=T)
cat("Click on Graphics Window to See Next Figure")
plot(density(mu2),type="l", lwd=2,col="darkgreen", 
	main="Marginal Posterior of Mu",xlab="Mu",ask="T")

lines(density(Mu[501:nsim]),lty=2, lwd=2, col="blue4")
legend(cex=.8,"topright",col=c("darkgreen",NA,"blue4"),lty=c(1,NA,2),lwd=2,
	legend=c("Direct Sampling from","Marginal Posterior","Gibbs Sampler"))


#######################################
# Draw sigma^2 directly from marginal #
#    posterior assuming non-          #
#      informative priors             #
#######################################
tau2<-rgamma(3500,(n-1)/2,(n-1)*var(y)/2)
sigma2_2<-1/tau2
mean(sigma2_2)
sd(sigma2_2)

par(ask=T)
cat("Click on Graphics Window to See Next Figure")
plot(density(sigma2_2),type="l", lwd=2,col="darkgreen", 
     main="Marginal Posterior of Sigma^2",xlab="Sigma^2",ask="T")

lines(density(Sigma[501:nsim]),lty=2, lwd=2, col="blue4")
legend(cex=.8,"topright",col=c("darkgreen",NA,"blue4"),lty=c(1,NA,2),lwd=2,
       legend=c("Direct Sampling from","Marginal Posterior","Gibbs Sampler"))

