
set.seed(100)
data = c(rnorm(1000,mean=2,sd=1),rnorm(400,mean=6,sd=1.5))

# linear
plot(sample(data), rep(0,length(data)), pch=4)

# with hist
hist(data, freq=FALSE, breaks=100)
points(sample(data), rep(0,length(data)), pch=4)

# with density
plot(density(data), col='red')
points(sample(data), rep(0,length(data)), pch=4)

# both
hist(data, freq=FALSE, breaks=100)
points(sample(data), rep(0,length(data)), pch=4)
lines(density(data), col='red')

normplot <- function(data,mu,sigma,lambda,color,...){
    x = seq(min(data),max(data),length=1000)
    y = lambda*dnorm(x,mu,sigma)
    lines(x,y,col=color,lwd=2)
}


#################
## With mclust ##
#################

library(mclust)
mixmodel = Mclust(data, G=2, model="V")  # G=2 components, "V" means unequal variances

mu = mixmodel$par$mean  ## should be approx [0,10]
sigma = sqrt(mixmodel$par$variance$sigmasq)  ## should be [1,2]

mclust1Dplot(data, parameters=mixmodel$parameters, z=mixmodel$z, what="density", col='red', new=TRUE, lwd=2, ylim=c(0,0.3))
hist(data, freq=FALSE, breaks=100, add=TRUE)
abline(v=mu, col='blue', lwd=2)
abline(v=c(mu+sigma,mu-sigma), col='blue', lty=2)


###################
## With mixtools ##
###################

# http://cran.r-project.org/web/packages/mixtools/mixtools.pdf

library(mixtools)
mixmodel = normalmixEM(data, k=2)  # k=2 composantes

mu = mixmodel$mu
sigma = mixmodel$sigma
post = mixmodel$posterior
lambda = mixmodel$lambda

# Default plot
plot(mixmodel, which=2)

# Custom plot
plot(density(data), col='red', ylim=c(0,0.3), xlim=c(-1.5,10), lwd=2)
normplot(data,mu[1],sigma[1],lambda[1],color='blue')
normplot(data,mu[2],sigma[2],lambda[2],color='green')
#points(sample(data), rep(0,length(data)), pch=4)
#abline(v=mu, col=c('blue','green'))
#abline(v=c(mu+sigma,mu-sigma), lty=2, col=c('blue','green'))
#lines(density(data), col='red', ylim=c(0,0.3), lwd=2)

plot(density(data), col='red', ylim=c(0,0.3), lwd=2)
x = seq(min(data),max(data),length=1000)
y1 = lambda[1]*dnorm(x,mu[1],sigma[1])
y2 = lambda[2]*dnorm(x,mu[2],sigma[2])
lines(x,y1,col='blue',lwd=2)
lines(x,y2,col='green',lwd=2)


##################
## Multivariate ##
##################

set.seed(1)
pop1 = mvrnorm(n=1000, c(2,2), matrix(c(1,0,0,1),2) )
pop2 = mvrnorm(n=1000, c(7,7), matrix(c(1.5,0,0,9),2) )
pop3 = mvrnorm(n=1000, c(1,8), matrix(c(10,0,0,5),2) )
data = rbind(pop1,pop2,pop3)

mvmixmodel = mvnormalmixEM(data, k=3)

par(mfrow=c(1,2))
plot(data,pch=4)
plot(mvmixmodel, which=2, pch=4)
