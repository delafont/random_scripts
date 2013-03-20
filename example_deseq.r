
# From DESeq's source code (tested)
estimate_size_factors <- function(counts) {
    loggeomeans = rowMeans(log(counts))
    apply( counts, 2, function(cnts){
        print(exp(log(cnts)-loggeomeans))
        exp( median( (log(cnts)-loggeomeans)[ is.finite(loggeomeans) ]))}
    )
}

# (tested)
estimate_mean <- function(counts){
    sf = estimate_size_factors(counts)
    qhat = t(t(counts)/sf)
    qpooled = counts$A/sf[1] + counts$B/sf[2]
    muA = sf[1]*qpooled
    muB = sf[2]*qpooled
    mu = data.frame(A=muA,B=muB,row.names=c('g1','g2','g3'))
}

# Neg bin example NB(16,10):
x = seq(100000)/100000
hist(rnbinom(x,mu=16,size=10),breaks=100)
abline(v=16,col='red')

counts = data.frame(A=c(4,12,9),B=c(4,3,4), row.names=c('g1','g2','g3'))
# Test if gene1 has same expression between conds A and B
mu = estimate_mean(counts)
sigma = mu # say var = mean
#f1 <- function(a){pnbinom(a,mu=mu[1,1], size=sigma[1,1])}
#f2 <- function(b){pnbinom(b,mu=mu[1,2], size=sigma[1,2])}
f1 <- function(a){ppois(a,1/mu[1,1])}
f2 <- function(b){ppois(b,1/mu[1,2])}
f <- function(a,b){f1(a)*f2(b)}
a = c()
for (i in 0:8) {a = c(a,f(i,8-i))}
sum(a[which(a<=a[5])]) / sum(a) # p-value for a=4, b=4

supercounts = 10*counts
# Test if gene1 has same expression between conds A and B
supermu = estimate_mean(supercounts)
supersigma = supermu # say var = mean
#g1 <- function(a){pnbinom(a,mu=supermu[1,1], size=supersigma[1,1])}
#g2 <- function(b){pnbinom(b,mu=supermu[1,2], size=supersigma[1,2])}
g1 <- function(a){ppois(a,1/supermu[1,1])}
g2 <- function(b){ppois(b,1/supermu[1,2])}
g <- function(a,b){g1(a)*g2(b)}
b = c()
for (i in 0:80) {b = c(b,g(i,80-i))}
sum(b[which(b<=b[41])]) / sum(b) # p-value for a=40, b=40

