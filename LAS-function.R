###----------------------------------------------###
###   R code for log-adjusted shrinkage prior    ###
###----------------------------------------------###

## This R code implements log-adjusted shrinkage (LAS) 
## prior for estimating multiple normal means.


## install packages 
library(GIGrvg)
library(MCMCpack)






###   LAS prior   ###
## INPUT
# y: vector of observations
# a: shape parameter for sparsity (default value is 1/n) 
# gam: shape parameter for log-adjustment (default value is 1/n) 
# L: order of iterated logarithm 
# mc: number of MCMC iterations
# burn: length of burn-in period 
## OUTPUT
# est: posterior means of signals
# theta: posterior samples of signals
# tau: posterior samples of tau (global scale parameter)

LAS <- function(y, a=NULL, gam=NULL, L=1, mc=3000, burn=1000, print=F){
  ## iterative log-function
  IL <- function(z, tt){ 
    val <- z
    if(tt>0){ 
      for(i in 1:tt){ val <- log(1+val) } 
    }
    return(val)
  }
  
  ## settings
  b <- 0
  n <- length(y)   # number of observations
  Delta <- 1/n    # hyperparameter of tau
  ep <- 10^(-10)
  if(is.null(a)){ a <- 1/n } 
  if(is.null(gam)){ gam <- 1/n } 
  
  ## MCMC box
  Th.pos <- matrix(NA, mc, n)
  Tau.pos <- rep(NA, mc)
  
  ## initial value
  theta <- y
  tau <- 1   # random effect variance
  U <- TT <- W <- rep(1, n)
  
  ##  MCMC iteration
  for(k in 1:mc){
    # Theta
    Kappa <- 1/(1+tau*U)
    mean_theta <- (1-Kappa)*y
    sd_theta <- sqrt(1-Kappa)
    theta <- rnorm(n, mean = mean_theta, sd = sd_theta)
    Th.pos[k,] <- theta
    # U
    a1 <- a - 1/2
    a2 <- theta^2/tau
    a3 <- 2*(TT+W)
    if(a1==0){ 
      a3[a3<ep] <- ep
      a2[a2<ep] <- ep 
    }
    if(a1<0){  a2[a2<ep] <- ep  }
    for(i in 1:n){  U[i] <- rgig(1, a1, a2[i], a3[i])  }
    # Tau
    vv <- rinvgamma(1, 1, 1/tau + 1/Delta^2 )   # latent variable for Half-Cauchy
    tau <- rinvgamma(1, n/2 + 1/2, 0.5*sum(theta^2/U) + 1/vv )
    Tau.pos[k] <- tau
    # TT
    TT <- rgamma(n, 1 + gam, 1+IL(U, L) )
    if(L>1){
      for(l in 1:L){
        TT <- rgamma(n, 1 + TT, 1+IL(U, L-l) )
      }
    }
    # W
    W <- rgamma(n, a + b, 1 + U)
    
    # print
    if(print){ if(round(k/100)==(k/100)) { print(k) } }
  }
  
  ## Summary
  sub <- (burn+1):mc
  Th.pos <- Th.pos[sub,]
  Tau.pos <- Tau.pos[sub]
  PM <- apply(Th.pos, 2, mean)
  Result <- list(est=PM, theta=Th.pos, tau=Tau.pos)
  return(Result)
}







###   adaptive LAS prior   ###
###  (gamma is estimated & L=1)  ###

## INPUT
# y: vector of observations
# a: shape parameter for sparsity (default value is 1/n) 
# mc: number of MCMC iterations
# burn: length of burn-in period 
## OUTPUT
# est: posterior means of signals
# theta: posterior samples of signals
# tau: posterior samples of tau (global scale parameter)
# gam: posteripr samples of gamma

aLAS <- function(y, a=NULL, mc=3000, burn=1000, K=50, print=F){
  ## functions
  # integrand of Riemann sum
  gfunc <- function(vx, db, dg){
    vy <- ( ( 1-exp(-vx) ) ^ (1-db) ) * ( ( 1+vx ) ^ (1+dg) ) 
    vy <- 1 / vy
    return(vy)
  }
  
  # Compute the upper/lower bounds of normalizing constant
  fbounds <- function(db, dg, k){
    N <- k^3	 
    U <- L <- 0
    U <- exp(1/k) * ( ( 1-exp(-1/k) )^(db) ) / db + ( ( 1-exp(-k) )^(db-1) ) / ( dg*((1+k)^(dg)) ) 
    L <- (1+1/k)^(-dg) * ( ( 1-exp(-1/k) )^(db) ) / db + (1+k)^(-dg) / dg 
    vx <- ( 1 + range(0,N-1) * (k^2-1) / N ) / k	
    U <- U + sum(gfunc(vx,db,dg)) * (k^2-1) / (k*N)
    vx <- ( 1 + range(1,N) * (k^2-1) / N ) / k	
    L <- L + sum(gfunc(vx,db,dg)) * (k^2-1) / (k*N)
    return(c(L, U))
  }
  
  # Sampling gamma by MH step with squeezing
  SampGam <-  function(da1, db1, dg0, a, dk0, n, ninc){ 
    dk <- dk0
    dg <- rgamma(1, da1, db1)
    fmis <- 1
    while(fmis){
      vbounds0 <- fbounds(a,dg0,dk)
      dwu <- dwl <-	n*a*( log(dg0)-log(dg) ) 
      dwu <- dwu + n*log( vbounds0[2] )
      dwl <- dwl + n*log( vbounds0[1] )
      vbounds <- fbounds(a,dg,dk)
      dwu <- dwu - n*log( vbounds[1] )
      dwl <- dwl - n*log( vbounds[2] )
      dwu <- exp(dwu)
      dwl <- exp(dwl)
      if(dwu > 1){ dwu <- 1 }
      if(dwl > 1){ dwl <- 1 }
      du <- runif(1)
      if(du < dwl){dacc <- 1; drej <- 0; fmis <- 0}                   # Accept
      else if(du > dwu){dg <- dg0; dacc <- 0; drej <- 1; fmis <- 0}    # Reject
      else{ dk <- dk + ninc}                                      # Missing (redo)  
    } 
    return(c(dg, dacc, drej))
  }
  
  ## settings
  b <- 0
  n <- length(y)   # number of observations
  Delta <- 1/n    # hyperparameter of tau
  ep <- 10^(-10)
  if(is.null(a)){ a <- 1/n }
  
  ## MCMC box
  Th.pos <- matrix(NA, mc, n)
  Tau.pos <- rep(NA, mc)
  Gam.pos <- rep(NA, mc)
  
  # initial value
  theta <- y
  tau <- 1   # random effect variance
  gam <- 0.5
  U <- TT <- W <- rep(1, n)
  
  # MCMC
  for(k in 1:mc){
    # Theta
    Kappa <- 1/(1+tau*U)
    mean_theta <- (1-Kappa)*y
    sd_theta <- sqrt(1-Kappa)
    theta <- rnorm(n, mean = mean_theta, sd = sd_theta)
    Th.pos[k,] <- theta
    # U
    a1 <- a - 1/2
    a2 <- theta^2/tau
    a3 <- 2*W
    if(a1==0){ 
      a3[a3<ep] <- ep
      a2[a2<ep] <- ep 
    }
    if(a1<0){  a2[a2<ep] <- ep  }
    for(i in 1:n){  U[i] <- rgig(1, a1, a2[i], a3[i])  }
    # Tau
    vv <- rinvgamma(1, 1, 1/tau + 1/Delta^2 )   # latent variable for HC
    tau <- rinvgamma(1, n/2 + 1/2, 0.5*sum(theta^2/U) + 1/vv )
    Tau.pos[k] <- tau
    # V
    V <- rgamma(n, 1 + gam, 1 + log(1+U) )
    # W
    W <- rgamma(n, V + a + b, 1 + U)
    # gam
    res <- SampGam(0.5+n*a, 0.5+sum(log(1+log(1+U))), gam, a, n, K, 10)
    gam <- min(res[1], 1)
    Gam.pos[k] <- gam
    
    # print
    if(print){ if(round(k/1000)==(k/1000)) { print(k) } }
  }
  
  ## Summary
  sub <- (burn+1):mc
  Th.pos <- Th.pos[sub,]
  Tau.pos <- Tau.pos[sub]
  Gam.pos <- Gam.pos[sub]
  PM <- apply(Th.pos, 2, mean)
  Result <- list(est=PM, theta=Th.pos, tau=Tau.pos, gam=Gam.pos)
  return(Result)
}












###  Three-parameter-Beta prior  ###
###     (default: HS prior)      ###
TPB <- function(y, a=1/2, b=1/2, mc=3000, burn=1000){
  ## settings
  p <- length(y)   # number of observations
  Delta <- 1/p    # hyperparameter of tau
  ep <- 10^(-10)
  
  ## initial value
  theta <- y
  tau <- 1   # random effect variance
  U <- W <- rep(1, p)
  
  ## MCMC box
  Theta.pos <- matrix(NA, mc, p)
  Tau.pos <- rep(NA, mc)
  
  ## MCMC
  for(k in 1:mc){
    # Theta
    Kappa <- 1/(1+tau*U)
    mean_theta <- (1-Kappa)*y
    sd_theta <- sqrt(1-Kappa)
    theta <- rnorm(p, mean = mean_theta, sd = sd_theta)
    Theta.pos[k,] <- theta
    # U
    a1 <- a - 1/2 
    a2 <- theta^2/tau
    a3 <- 2*W
    if(a1==0){ 
      a3[a3<ep] <- ep
      a2[a2<ep] <- ep 
    }
    if(a1<0){ a2[a2<ep] <- ep }
    for(i in 1:p){ U[i] <- rgig(1, a1, a2[i], a3[i]) }
    
    # Tau
    vv <- rinvgamma(1, 1, 1/tau + 1/Delta^2 )   # latent variable for HC
    tau <- rinvgamma(1, p/2 + 1/2, 0.5*sum(theta^2/U) + 1/vv )
    Tau.pos[k] <- tau
    # W
    W <- rgamma(p, a + b, 1 + U)
  }
  
  ## Summary
  omit <- 1:burn
  Theta.pos <- Theta.pos[-omit,]
  Tau.pos <- Tau.pos[-omit]
  PM <- apply(Theta.pos, 2, mean)
  Result <- list(est=PM, theta=Theta.pos, tau=Tau.pos)
  return(Result)
}







