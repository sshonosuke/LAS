## settings
rm(list=ls())
set.seed(1)
mc <- 10000    # MCMC length
bn <- 5000    # length of burn-in 

## load functions and packages
source("LAS-function.R")
library(NormalBetaPrime)


## load dataset
load("prostz.Rdata")

## data 
Y <- prostz
n <- length(Y)


## log-adjusted shrinkage prior (LAS)
LS <- LAS(Y, a=1/n^2, gam=1/n, mc=mc, burn=bn) 

## adaptive LAS (aLAS)
aLS <- aLAS(Y, a=1/n, mc=mc, burn=bn)

## iteratively log-adjusted shrinkage prior (ILAS)
ILS <- LAS(Y, a=1/n, gam=1/n, L=3, mc=mc, burn=bn)

## horseshoe (HS)
HS <- TPB(Y, a=1/2, b=1/2, mc=mc, burn=bn) 
 
## horseshoe plus (HS+)
HSP <- hsplus.normalmeans(Y, tau.est="truncatedCauchy", max.steps=mc, burnin=bn,  var.select="threshold")

## normal prime beta (NPB)
NPB <- nbp.normalmeans(Y, a.est="fixed", max.steps=mc, burnin=bn,  var.select="threshold")

## Dirichlet Laplace (DL)
DL <- dl.normalmeans(Y, tau.est="truncatedCauchy", max.steps=mc, burnin=bn,  var.select="threshold")


## point estimate
Lab <- c("Z", "LAS", "aLAS", "ILAS", "HS", "HS+", "NPB", "DL")
Est <- cbind(Y, LS$est, aLS$est, ILS$est, HS$est, HSP$theta.hat, NPB$theta.hat, DL$theta.hat)
dimnames(Est)[[2]] <- Lab
head(Est)



## posterior predictive loss (PPL)
PPL <- function(pos){
  pos2 <- pos + rnorm(prod(dim(pos)))
  v1 <- mean( apply(pos2, 2, var) ) 
  v2 <- mean( (Y - apply(pos2, 2, mean))^2 )*n/(n+1)
  return( c(v1, v2) )
}

Loss <- matrix(NA, 7, 3)
Loss[1,-3] <- PPL(LS$theta)
Loss[2,-3] <- PPL(aLS$theta)
Loss[3,-3] <- PPL(ILS$theta)
Loss[4,-3] <- PPL(HS$theta)
Loss[5,-3] <- c(1 + mean(HSP$theta.var), mean( (Y - HSP$theta.hat)^2 )*n/(n+1))
Loss[6,-3] <- c(1 + mean(NPB$theta.var), mean( (Y - NPB$theta.hat)^2 )*n/(n+1))
Loss[7,-3] <- c(1 + mean(DL$theta.var), mean( (Y - DL$theta.hat)^2 )*n/(n+1))
Loss[,3] <- apply(Loss[,1:2], 1, sum)

dimnames(Loss)[[1]] <- Lab
dimnames(Loss)[[2]] <- c("Variance", "Squared error", "PPL")
Loss



## scatter plot of posterior means
plot(Y, LS$est, pch=16, cex=1, xlab="Observation", ylab="Posterior mean")
abline(0, 1)
abline(h=0, lty=2)
points(Y, NPB$theta.hat, col="red", pch=16, cex=1)
points(Y, HS$est, col="green", pch=16, cex=1)
points(Y, HSP$theta.hat, col="blue", pch=16, cex=1)
legend("topleft", legend=c("LAS", "NPB", "HS", "HS+"), col=c("black", "red", "green", "blue"), pch=16)



## trace plot of MCMC samples
mtheta1 <- apply(LS$theta, 1, mean)
mtheta2 <- apply(aLS$theta, 1, mean)
mtheta3 <- apply(ILS$theta, 1, mean)
mtheta4 <- apply(HS$theta, 1, mean)

par(mfcol=c(2,2))
plot(mtheta1, type="l", xlab="iteration", ylab="", main="LAS")
plot(mtheta2, type="l", xlab="iteration", ylab="", main="aLAS")
plot(mtheta3, type="l", xlab="iteration", ylab="", main="ILAS")
plot(mtheta4, type="l", xlab="iteration", ylab="", main="HS")






