## ----eval=FALSE---------------------------------------------------------------
#  sample_diploid <- function(n, lambda, alpha, p){
#    reads <- vector("list", length = n)
#    geno <- sample(c(2,1,0),size = n,replace = TRUE, prob = c((1-p)^2,2*(1-p)*p,p^2))
#    N <- rpois(n,lambda)
#    for(i in 1:n){
#      nR <- geno[i]
#      pR <- nR/2+(1-nR)*alpha
#      reads[[i]] <- sample(c(1,0),size = N[i],replace = TRUE,prob = c(pR,1-pR))
#    }
#    return(reads)
#  }

## ----eval=FALSE---------------------------------------------------------------
#  seqem_diploid <- function(reads, iters, epsilon){
#    alpha0 <- 0.01
#    p0 <- 0.2
#    n <- length(reads)
#    N <- numeric(n)
#    X <- numeric(n)
#    for(i in 1:n){
#      N[i] <- length(reads[[i]])
#      X[i] <- sum(reads[[i]]==0)
#    }
#  
#    for(j in 1:iters){
#      f0 <- f1 <- f2 <- numeric(n)
#  
#      for(i in 1:n){
#        f0[i] <- (1-alpha0)^X[i]*alpha0^(N[i]-X[i])*p0^2
#        f1[i] <- (1/2)^N[i]*2*p0*(1-p0)
#        f2[i] <- alpha0^X[i]*(1-alpha0)^(N[i]-X[i])*(1-p0)^2
#      }
#      f <- f0+f1+f2
#      f0 <- f0/f
#      f1 <- f1/f
#      f2 <- f2/f
#  
#      X_VV <- sum(X*f0)
#      X_RR <- sum(X*f2)
#      N_VV <- sum(N*f0)
#      N_RR <- sum(N*f2)
#      S_VV <- sum(f0)
#      S_RV <- sum(f1)
#  
#      alpha1 <- (N_VV-X_VV+X_RR)/(N_VV+N_RR)
#      p1 <- (2*S_VV+S_RV)/(2*n)
#  
#      if((alpha1<0.5) && (abs(alpha1-alpha0)<epsilon) && (abs(p1-p0)<epsilon)){
#        return(c(alpha1,p1))
#      }
#      else{
#        alpha0 <- alpha1
#        p0 <- p1
#      }
#    }
#    return(c(-1,-1))
#  }

## -----------------------------------------------------------------------------
trials_diploid <- function(m,n,lambda,alpha,p){
  es_alpha <- numeric(m)
  es_p <- numeric(m) 
  index <- numeric(m) 
  
  for(i in 1:m){
    sample <- sample_diploid(n,lambda,alpha,p)
    em <- seqem_diploid(sample,200,1e-8)
    es_alpha[i] <- em[1]
    es_p[i] <- em[2]
    if(em[1]>0.1 || em[1]==-1) index[i] <- FALSE 
    else index[i] <- TRUE
  }
  es <- rbind(es_alpha,es_p)
  es <- es[,index==TRUE]
  return(es)
}

## ----eval=TRUE----------------------------------------------------------------
library(StatComp21017)
tr <- trials_diploid(1000,50,10,0.01,0.2)
apply(tr,1,mean)
apply(tr,1,sd)

## -----------------------------------------------------------------------------
sample_tetraploid <- function(n,lambda,alpha,p){
  reads <- vector("list",length = n)
  geno <- sample(c(4,3,2,1,0),size=n,replace=T,prob=c((1-p)^4,4*(1-p)^3*p,6*(1-p)^2*p^2,4*(1-p)*p^3,p^4))
  N <- rpois(n,lambda)
  for(i in 1:n){
    nR <- geno[i]
    pR <- nR/4+(1-nR/2)*alpha
    reads[[i]] <- sample(c(1,0),size=N[i],replace=T,prob=c(pR,1-pR))
  }
  return(reads)
}

## -----------------------------------------------------------------------------
seqem_tetraploid <- function(reads, iters, epsilon){
  alpha0 <- 0.01
  p0 <- 0.2 
  n <- length(reads)
  N <- numeric(n)
  X <- numeric(n)
  for(i in 1:n){
    N[i] <- length(reads[[i]])
    X[i] <- sum(reads[[i]]==0)
  }
  
  for(j in 1:iters){
    f0 <- f1 <- f2 <- f3 <- f4 <- numeric(n)
    
    for(i in 1:n){
      f0[i] <- (1-alpha0)^X[i]*alpha0^(N[i]-X[i])*p0^4
      f1[i] <- (3/4-alpha0/2)^X[i]*(1/4+alpha0/2)^(N[i]-X[i])*4*(1-p0)*p0^3
      f2[i] <- (1/2)^N[i]*6*(1-p0)^2*p0^2
      f3[i] <- (1/4+alpha0/2)^X[i]*(3/4-alpha0/2)^(N[i]-X[i])*4*(1-p0)^3*p0
      f4[i] <- alpha0^X[i]*(1-alpha0)^(N[i]-X[i])*(1-p0)^4
    }
    f <- f0+f1+f2+f3+f4
    f0 <- f0/f 
    f1 <- f1/f 
    f2 <- f2/f 
    f3 <- f3/f 
    f4 <- f4/f 
    
    X0 <- sum(X*f0);X1 <- sum(X*f1);X3 <- sum(X*f3);X4 <- sum(X*f4)
    N0 <- sum(N*f0);N1 <- sum(N*f1);N3 <- sum(N*f3);N4 <- sum(N*f4)
    S0 <- sum(f0);S1 <- sum(f1);S2 <- sum(f2);S3 <- sum(f3)
    
    alpha1 <- uniroot(function(x) (X4+N0-X0)/x+(X0+N4-X4)/(x-1)+(X3+N1-X1)/(x+1/2)+(X1+N3-X3)/(x-3/2),lower = 0,upper = 0.1,extendInt = "yes")$root
    p1 <- (4*S0+3*S1+2*S2+S3)/(4*n)
    
    if((abs(alpha1-alpha0)<epsilon) && (abs(p1-p0)<epsilon)){
      return(c(alpha1,p1))
    }
    else{
      alpha0 <- alpha1
      p0 <- p1 
    }
  }
  return(c(-1,-1)) 
}

## -----------------------------------------------------------------------------
trials_tetraploid <- function(m,n,lambda,alpha,p){
  es_alpha <- numeric(m)
  es_p <- numeric(m) 
  index <- numeric(m) 
  
  for(i in 1:m){
    sample <- sample_tetraploid(n,lambda,alpha,p)
    em <- seqem_tetraploid(sample,300,1e-8)
    es_alpha[i] <- em[1]
    es_p[i] <- em[2]
    if(em[1]>0.1 || em[1]==-1) index[i] <- FALSE 
    else index[i] <- TRUE
  }
  es <- rbind(es_alpha,es_p)
  es <- es[,index==TRUE]
  return(es)
}

## ----eval=TRUE----------------------------------------------------------------
tr <- trials_tetraploid(1000,100,8,0.01,0.1)
apply(tr,1,mean)
apply(tr,1,sd)

