#################
## With mclust ##
#################

library(mclust)
data = c(rnorm(1000),rnorm(400,sd=2,mean=10))
clust = Mclust(data,G=2,model="V") #G=2 populations, "V" means unequal variances

mu = clust$par$mean  ## should be approx [0,10]
sigma = sqrt(clust$par$variance$sigmasq)  ## should be [1,2]

mclust1Dplot(data,parameters=clust$parameters,z=clust$z,what="density",col='red',new=TRUE)
hist(data,freq=FALSE,breaks=100,add=TRUE)
abline(v=mu,col='blue')
abline(v=c(mu+sigma,mu-sigma),col='blue',lty=2)

###################
## With mixtools ##
###################

library(mixtools)
mixmdl = normalmixEM(data)
plot(mixmdl,which=2)
lines(density(data), lty=2, lwd=2)

