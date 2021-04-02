rm(list=ls())
set.seed(1)


# load functions and packages
source("LAS-function.R")
library(NormalBetaPrime)

quant <- function(x){ quantile(x, prob=c(0.025, 0.975)) }


# scenario selection 
sc <- 1    # 1, 2 or 3
C <- 9    # 6 or 9 (signal location) 
om <- 0.2   # 0.1, 0.2 or 0.3


# settings
n <- 200    # number of observations
nn <- round(n*om)



# generation of true signals 
th <- rep(0, n)

if(sc==1){ 
  th[1:(nn/2)] <- C
  th[(nn/2+1):nn] <- (-C/2)
}

if(sc==2){ 
  th[1:(nn/2)] <- rnorm(nn/2, C, 1)
  th[(nn/2+1):nn] <- rnorm(nn/2, -C/2, 1)
}

if(sc==3){ 
  th[1:nn] <- rnorm(nn, 0, sqrt(2*log(n)))
}

# generation of observations
Y <- th + rnorm(n)


## estimation 
Lab <- c("LAS", "aLAS", "ILAS", "HS", "HS+", "NPB", "DL")
L <- length(Lab)

# LAS prior 
LS <- LAS(Y, a=1/n, gam=1/n, mc=2000, burn=1000)  

# adaptive LAS prior 
aLS <- aLAS(Y, a=1/n, mc=2000, burn=1000)

# iterated LAS prior
ILS <- LAS(Y, a=1/n, gam=1/n, L=3, mc=2000, burn=1000)  

# Horseshoe prior
HS <- TPB(Y, a=1/2, b=1/2, mc=2000)   

# Horseshoe+ prior
HSP <- hsplus.normalmeans(Y, tau.est="truncatedCauchy", max.steps=2000, burnin=1000,  var.select="threshold")

# Normal-prime-beta prior
NPB <- nbp.normalmeans(Y, a.est="fixed", max.steps=2000, burnin=1000,  var.select="threshold")

# Dirichlet-Laplace prior
DL <- dl.normalmeans(Y, tau.est="truncatedCauchy", max.steps=2000, burnin=1000,  var.select="threshold")



# L2-loss
Est <- cbind(LS$est, aLS$est, ILS$est, HS$est, HSP$theta.hat, NPB$theta.hat, DL$theta.hat)
loss <- apply((Est-th)^2, 2, sum) 
names(loss) <- Lab
loss


# Interval 
CI <- list(apply(LS$theta, 2, quant), apply(aLS$theta, 2, quant), apply(ILS$theta, 2, quant), 
           apply(HS$theta, 2, quant), HSP$theta.intervals, NPB$theta.intervals, DL$theta.intervals)
CP <- c()
AL <- c()
IS <- c()
for(l in 1:L){
  CP[l] <- mean( CI[[l]][1,]<th & CI[[l]][2,]>th )
  AL[l] <- mean( CI[[l]][2,]-CI[[l]][1,] )
  IS[l] <- AL[l] + mean( 40*ifelse(th>CI[[l]][2,], 1, 0)*(th-CI[[l]][2,]) + 40*ifelse(th<CI[[l]][1,], 1, 0)*(CI[[l]][1,]-th) )
}

names(CP) <- names(AL) <- names(IS) <- Lab
CP*100
AL
IS
