## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
library(knitr)
library(xtable)

## -----------------------------------------------------------------------------
lm <- lm(Petal.Length~Petal.Width,data = iris)
plot(lm)

## -----------------------------------------------------------------------------
irisset <- subset(iris,Species=="setosa")
boxplot(irisset[,1:4],main="setosa",col = c(1,2,3,4))

## -----------------------------------------------------------------------------
knitr::kable(head(iris[,1:4]))

## ----xtable,results='asis'----------------------------------------------------
print(xtable(head(iris[,1:4])),type = "html")

## -----------------------------------------------------------------------------
set.seed(1234)

# generate n random number from Rayleigh distribution, the parameter is sigma
rRayleigh <- function(n,sigma){
  u <- runif(n)
  x <- sqrt(-2*sigma^2*log(1-u))
  return(x)
}

# generate sample with sigma=5
n <- 1000
sigma <- 5
x <- rRayleigh(n,sigma)

# check generation algorithm is correct
hist(x,probability = T,main = expression(sigma==5))
y <- seq(0,20,.01)
lines(y,y/25*exp(-y^2/50))

# check the mode of samples is close to the sigma
hist(x=rRayleigh(n,sigma=2),main = expression(sigma==2))
hist(x=rRayleigh(n,sigma=4),main = expression(sigma==4))
hist(x=rRayleigh(n,sigma=5),main = expression(sigma==5))
hist(x=rRayleigh(n,sigma=10),main = expression(sigma==10))

## -----------------------------------------------------------------------------
set.seed(1234)
n <- 1000  # sample size = 1000

# generate x1,x2,u
x1 <- rnorm(n,mean = 0,sd=1)
x2 <- rnorm(n,mean = 3,sd=1)
u <-runif(n)

# generate x with mixture distribution
rmixture <- function(p1=0.5,p2){
  x <- vector(length = n)
  for(i in 1:n){
    if(u[i]<=p1) x[i] <- x1[i]
    if(u[i]>p1) x[i] <- x2[i]
  }
  return(x)
}

# generate samples with p1=0.5
x <- rmixture(0.75,0.25)

# hist plot with density superimposed
hist(x,probability = T,main="p1=0.75",ylim = c(0,0.3))
lines(density(x))

# vary p1
for(p1 in seq(0.1,0.9,0.1)){
  x <- rmixture(p1,1-p1)
  hist(x,probability = T,main=c("p1=",p1),ylim = c(0,0.4))
  lines(density(x))
}

## -----------------------------------------------------------------------------
n <- 1000
rPG <- function(lambda,t,alpha,beta){
  x <- vector(length = n)
  for(i in 1:n){
    Nt <- rpois(1,lambda=lambda*t) # Nt ~ Poi(lambda*t)
    y <- rgamma(Nt,alpha,beta)
    x[i] <- sum(y)
  }
  return(x)
}

# compare the estimated mean and variance with the theoretical values
results <- matrix(0,nrow = 4,ncol = 5)
xx <- matrix(0,nrow = 4, ncol = n) # original samples
colnames(results) <- c("parameter","mean","true mean","var","true var")

# several choices of the parameters
Lambda <- c(1,1,3,5)
Alpha <- c(1,1,2,3)
Beta <- c(1,2,1,0.5)
for(i in 1:4){
  xx[i,] <- rPG(Lambda[i],10,Alpha[i],Beta[i])
}

results[,2] <- round(rowMeans(xx),4)
results[,3] <- Lambda*10*Alpha/Beta
results[,5] <- Lambda*10*Alpha*(1+Alpha)/(Beta^2)
for(i in 1:4){
  results[i,1] <- paste("lambda=",Lambda[i],",alpha=",Alpha[i],",beta=",Beta[i])
  results[i,4] <- round(var(xx[i,]),4)
}

library(knitr)
knitr::kable(results)

## -----------------------------------------------------------------------------
set.seed(121)
x <- seq(.1,.9,length=9) # for some x in (0,1)
m <- 10000
y <- runif(m)
cdf <- numeric(length(x)) # F(x)

for(i in 1:length(x)){
  if(x[i]<=0 || x[i] >=1) cdf[i]<-0
  else{
    g <- 30*x[i]^3*y^2*(rep(1,m)-x[i]*y)^2
    cdf[i] <- mean(g)
  }
}
result <- rbind(x,cdf,pbeta(x,3,3))
rownames(result) <- c('x','cdf',"pbeta")
knitr::kable(result)

## -----------------------------------------------------------------------------
# generate samples
MC.Rayleigh <- function(sigma, m=10000, antithetic=TRUE){
  u <- runif(m)
  if(!antithetic) v <- runif(m) # independent
  else v <- 1-u
  x <- (sqrt(-2*sigma^2*log(1-u))+sqrt(-2*sigma^2*log(1-v)))/2
  return(x)
}

sigma <- c(0.5,1,2,5,10,15)
var_re <- numeric(length(sigma)) # reduction of variance
for(i in 1:length(sigma)){
  MC1 <- MC.Rayleigh(sigma[i],m=1000,antithetic = FALSE)
  MC2 <- MC.Rayleigh(sigma[i],m=1000,antithetic = TRUE)
  var_re[i] <- (var(MC1)-var(MC2))/var(MC1)
}

result2 <- rbind(sigma,var_re)
knitr::kable(result2)

## -----------------------------------------------------------------------------
m <- 10000
theta.hat <- var <- numeric(2)

g <- function(x){
  x^2*exp(-x^2/2)*(x>1)/(sqrt(2*pi))
}


# f1
x <- rnorm(m,mean = 1.5)
g_f <- g(x)/dnorm(x,mean = 1.5)
theta.hat[1] <- mean(g_f)
var[1] <- var(g_f)

# f2
u <- runif(m)
x <- 1/(1-u) # inverse transformation
g_f <- g(x)*x^2
theta.hat[2] <- mean(g_f)
var[2] <- var(g_f)

theta.hat
var

## -----------------------------------------------------------------------------
integrate(g,1,upper = Inf) # true value
#plot
y <- seq(1,10,.01)
plot(y,g(y),type = "l",main = "",ylab = "",ylim = c(0,2),col="red")
lines(y,dnorm(y,mean = 1.5),lty=2,col=3)
lines(y,1/(y^2),lty=3,col=4)

## -----------------------------------------------------------------------------
set.seed(200)
m <- 10000
g <- function(x){
  x^2*exp(-x^2/2)*(x>1)/(sqrt(2*pi))
}


# f1
x <- rnorm(m,mean = 1.5)
g_f <- g(x)/dnorm(x,mean = 1.5)
estimate <- mean(g_f)
true <- integrate(g,1,upper = Inf)$value # true value
c(estimate,true)

## -----------------------------------------------------------------------------
# Symmetric t-intervals
n <- 20 # sample size
alpha <- 0.05
set.seed(200)
CL <- replicate(1000,expr = { 
  x <- rchisq(n,df=2)
  L <- mean(x)-qt(1-alpha/2,df=n-1)*sd(x)/sqrt(n)   # Left end of interval
  R <-  mean(x)+qt(1-alpha/2,df=n-1)*sd(x)/sqrt(n)  # Right end of interval
  UCL <- (n-1)*var(x)/qchisq(alpha,df=n-1) 
  c(L,R,UCL)
})

# Calculate the empirical confidence level
m <- 0
for(i in 1:1000){
  if(CL[1,i]<=2&&CL[2,i]>=2) m <- m+1
}
m/1000

mean(CL[3,]>4) # The variance range

## -----------------------------------------------------------------------------
# Type I error
m <- 10000 
n <- 20    # sample size
alpha <- 0.05
p_chi <- p_unif <- p_expr <- numeric(m) # p values for each trial repeated

# three distributions
set.seed(123)
for(i in 1:m){
  chi <- rchisq(n,df=1)
  unif <- runif(n,0,2)
  expr <- rexp(n ,rate = 1)
  test1 <- t.test(chi,alternative="two.sided",mu=1,conf.level=0.95)
  test2 <- t.test(unif,alternative="two.sided",mu=1,conf.level=0.95)
  test3 <- t.test(expr,alternative="two.sided",mu=1,conf.level=0.95)
  p_chi[i] <- test1$p.value
  p_unif[i] <- test2$p.value
  p_expr[i] <- test3$p.value
}
p_chi.hat <- mean(p_chi < alpha)
p_unif.hat <- mean(p_unif < alpha)
p_expr.hat <- mean(p_expr < alpha)
c(p_chi.hat,p_unif.hat,p_expr.hat)

## -----------------------------------------------------------------------------
mat <-
  matrix(c(6510, 3490, 10000, 6760, 3240, 10000, 13270, 6730, 20000), 3, 3,
         dimnames = list(
           c("Rejected", "Accepted", "total"),
           c("method A", "method B", "total")
         ))
mat

## ----eval=FALSE---------------------------------------------------------------
#  nn <- c(10,20,30,50,100,500)  # sample size
#  alpha <- 0.05                 # significance level
#  d <- 2                        # The dimension of a random variable
#  b0 <- qchisq(1-alpha,df=d*(d+1)*(d+2)/6)*6/nn  # Each sample size threshold vector
#  
#  # Calculate multivariate sample skewness statistics
#  mul.sk <- function(x){
#    n <- nrow(x)
#    xbar <- colMeans(x)
#    sigma.hat <- (n-1)/n*cov(x) # MLE
#    V <- solve(sigma.hat)
#  
#    b <- 0
#    for(i in 1:nrow(x)){
#      for(j in 1:nrow(x)){
#        b <- b+((x[i,]-xbar)%*%V%*%(x[j,]-xbar))^3
#      }
#    }
#    return(b/(n^2))
#  }
#  
#  # calculate an empirical estimate of the first type of error
#  library(mvtnorm)
#  set.seed(200)
#  p.reject <- vector(mode = "numeric",length = length(nn))
#  
#  m <- 1000
#  
#  for(i in 1:length(nn)){
#    mul.sktests <- vector(mode = "numeric",length = m)
#    for(j in 1:m){
#      data <- rmvnorm(nn[i],mean = rep(0,d))
#      mul.sktests[j] <- as.integer(mul.sk(data)>b0[i])
#    }
#    p.reject[i] <- mean(mul.sktests)
#  }
#  p.reject

## ----eval=FALSE---------------------------------------------------------------
#  summ <- rbind(nn,p.reject)
#  rownames(summ) <- c("n","estimate")
#  knitr::kable(summ)

## ----eval=FALSE---------------------------------------------------------------
#  alpha <- 0.1
#  n <- 30      # sample size
#  m <- 2000
#  epsilon <- c(seq(0,0.15,0.01),seq(0.15,1,0.05))
#  N <- length(epsilon)
#  power <- vector(mode = "numeric",length = N)
#  b0 <- qchisq(1-alpha,df=d*(d+1)*(d+2)/6)*6/n  #Critical values
#  
#  # Separate power for this column of epsilon
#  for(j in 1:N){
#    e <- epsilon[j]
#    mul.sktests <- numeric(m)
#    for(i in 1:m){
#      # Generate mixed distribution
#      u <- sample(c(1,0),size = n,replace = T,prob = c(1-e,e))
#      data1 <- rmvnorm(n,sigma = diag(1,d))
#      data2 <- rmvnorm(n,sigma = diag(100,d))
#      data <- u*data1+(1-u)*data2
#      mul.sktests[i] <- as.integer(mul.sk(data)>b0)
#    }
#    power[j] <- mean(mul.sktests)
#  }
#  
#  # Plot the power function
#  plot(epsilon,power,type="b",xlab=bquote(epsilon),ylim=c(0,1))
#  abline(h=0.1,lty=3,col="lightblue")
#  se <- sqrt(power*(1-power)/m)
#  lines(epsilon,power-se,lty=3)
#  lines(epsilon,power+se,lty=3)

## -----------------------------------------------------------------------------
library(boot);library(bootstrap)
data(scor)

set.seed(200)
lambda_hat <- eigen(cov(scor))$values
theta_hat <- lambda_hat[1]/sum(lambda_hat)
B <- 2000     
boot.theta <- function(x,index){
  lambda <- eigen(cov(x[index,]))$values
  theta <- lambda[1]/sum(lambda)
  return(theta)
}

bootstrap_result <- boot(data = scor,statistic = boot.theta,R=B)
theta_b <- bootstrap_result$t
bias_boot <- mean(theta_b)-theta_hat   # bias
se_boot <- sqrt(var(theta_b))          # se
c(bias_boot=bias_boot,se_boot=se_boot)

## -----------------------------------------------------------------------------
# Jackknife
n <- nrow(scor)
theta_jack <- rep(0,n)
for(i in 1:n){
  x <- scor[-i,]
  lambda <- eigen(cov(x))$values
  theta_jack[i] <- lambda[1]/sum(lambda)
}
bias_jack <- (n-1)*(mean(theta_jack)-theta_hat) # bias
se_jack <- (n-1)*sqrt(var(theta_jack)/n)        # se
c(bias_jack=bias_jack,se_jack=se_jack)

## -----------------------------------------------------------------------------
ci <- boot.ci(bootstrap_result,type = c("perc","bca"))
ci

## -----------------------------------------------------------------------------
set.seed(200);m <- 2000
boot.skewness <- function(x,i){
  a <- sum((x[i]-mean(x[i]))^3)/length(x[i])
  b <- (sum((x[i]-mean(x[i]))^2)/length(x[i]))^(3/2)
  return(a/b)
}
ci.norm <- ci.basic <- ci.perc <- matrix(0,m,2)

# normal population
skew_normal <- 0
for(i in 1:m){
  xx <- rnorm(n=100,0,1)
  boot_res <- boot(data=xx,statistic=boot.skewness,R=999) # bootstrap
  ci <- boot.ci(boot_res,type=c("norm","basic","perc"))
  ci.norm[i,] <- ci$norm[2:3]
  ci.basic[i,] <- ci$basic[4:5]
  ci.perc[i,] <- ci$percent[4:5]
}

# Confidence interval coverage rate
coverage_rate <- c(mean(ci.norm[,1]<=skew_normal & ci.norm[,2]>=skew_normal),mean(ci.basic[,1]<=skew_normal & ci.basic[,2]>=skew_normal),mean(ci.perc[,1]<=skew_normal & ci.perc[,2]>=skew_normal))

# Percentage of misses on the left
miss_left <- c(mean(ci.norm[,1]>skew_normal),mean(ci.basic[,1]>skew_normal),mean(ci.perc[,1]>skew_normal))

# Proportion of misses on the right
miss_right <- c(mean(ci.norm[,2]<skew_normal),mean(ci.basic[,2]<skew_normal),mean(ci.perc[,2]<skew_normal))

cat('For normal populations, the estimate of the coverage probabilities of the normal bootstrap CI =',coverage_rate[1],",the basic CI =",coverage_rate[2],"and the percentile CI =",coverage_rate[3])

normal <- cbind(coverage_rate,miss_left,miss_right)
rownames(normal) <- c("normal","basic","percentile")
knitr::kable(normal)

## -----------------------------------------------------------------------------
# chi^2(5) distributions 
ci.norm <- ci.basic <- ci.perc <- matrix(0,m,2)
skew_chi <- sqrt(8/5)
for(i in 1:m){
  xx <- rchisq(n=100,df=5)
  boot_res <- boot(data=xx,statistic=boot.skewness,R=999) # bootstrap
  ci <- boot.ci(boot_res,type=c("norm","basic","perc"))
  ci.norm[i,] <- ci$norm[2:3]
  ci.basic[i,] <- ci$basic[4:5]
  ci.perc[i,] <- ci$percent[4:5]
}

coverage_rate <- c(mean(ci.norm[,1]<=skew_chi & ci.norm[,2]>=skew_chi),mean(ci.basic[,1]<=skew_chi & ci.basic[,2]>=skew_chi),mean(ci.perc[,1]<=skew_chi & ci.perc[,2]>=skew_chi))

# Percentage of misses on the left
miss_left <- c(mean(ci.norm[,1]>skew_chi),mean(ci.basic[,1]>skew_chi),mean(ci.perc[,1]>skew_chi))

miss_right <- c(mean(ci.norm[,2]<skew_chi),mean(ci.basic[,2]<skew_chi),mean(ci.perc[,2]<skew_chi))

cat('For chi-square distributions, the estimate of the coverage probabilities of the normal bootstrap CI =',coverage_rate[1],",the basic CI =",coverage_rate[2],"and the percentile CI =",coverage_rate[3])

chi <- cbind(coverage_rate,miss_left,miss_right)
rownames(chi) <- c("normal","basic","percentile")
knitr::kable(chi)

## -----------------------------------------------------------------------------
# Generate two related samples
library(MASS)
set.seed(200)
data <- mvrnorm(n=25,mu=c(1,2),Sigma = matrix(c(2,1,1,3),nrow=2,ncol=2))
x <- data[,1]
y <- data[,2]

# Permutation test
R <- 2000-1
z <- c(x,y)  # pooled sample
n1 <- length(x);n2 <- length(y)
K <- 1:(n1+n2)
reps <- vector(mode = "numeric",length = R)
rho_0 <- cor(x,y,method = "spearman")

for(i in 1:R){
  k <- sample(K, size = n1,replace = F)
  x1 <- z[k]
  y1 <- z[-k]
  reps[i] <- cor(x1,y1,method = "spearman")
}
p_permu <- mean(c(rho_0,reps)>=rho_0)
p_cortest <- cor.test(x,y)$p.value
round(c(p_permu,p_cortest),4)

# independent sample
data <- mvrnorm(n=25,mu=c(1,0),Sigma = matrix(c(2,0,0,3),nrow = 2))
x <- data[,1]
y <- data[,2]

# permutation test 
z <- c(x,y)
reps <- vector(mode = "numeric",length = R)
rho_0 <- cor(x,y,method = "spearman")

for(i in 1:R){
  k <- sample(K, size = n1,replace = F)
  x1 <- z[k]
  y1 <- z[-k]
  reps[i] <- cor(x1,y1,method = "spearman")
}
p_permu <- mean(abs(c(rho_0,reps))>=abs(rho_0))
p_cortest <- cor.test(x,y)$p.value
round(c(p_permu,p_cortest),4)

## ----eval=FALSE---------------------------------------------------------------
#  library(RANN)
#  library(boot)
#  library(energy)
#  library(Ball)
#  library(MASS)
#  
#  # NN
#  Tn <- function(z, index, sizes,k) {
#    n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
#    if(is.vector(z)) z <- data.frame(z,0);
#    z <- z[index, ];
#    NN <- nn2(data=z, k=k+1)
#    block1 <- NN$nn.idx[1:n1,-1]
#    block2 <- NN$nn.idx[(n1+1):n,-1]
#    i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
#    (i1 + i2) / (k * n)
#  }
#  
#  eqdist.nn <- function(z,sizes,k){
#    boot.obj <- boot(data=z,statistic=Tn,R=R, sim = "permutation", sizes = sizes,k=k)
#    ts <- c(boot.obj$t0,boot.obj$t)
#    p.value <- mean(ts>=ts[1])
#    list(statistic=ts[1],p.value=p.value)
#  }

## ----eval=FALSE---------------------------------------------------------------
#  k <- 3;R <- 500-1;m <- 200
#  
#  # setting
#  n1 <- n2 <- 20
#  N <- c(n1,n2)
#  set.seed(123)
#  
#  # power
#  p.values <- matrix(NA,m,3)
#  for(i in 1:m){
#    x <- mvrnorm(n1,c(0,0),diag(1,2))
#    y <- mvrnorm(n2,c(0,0),diag(3,2))
#    z <- rbind(x,y)
#    p.values[i,1] <- eqdist.nn(z,N,k)$p.value
#    p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
#    p.values[i,3] <- bd.test(x,y,num.permutations=R,seed=i*123)$p.value
#  }
#  alpha <- 0.1;
#  pow <- colMeans(p.values<alpha)
#  pow

## ----eval=FALSE---------------------------------------------------------------
#  k <- 3;R <- 500-1;m <- 200
#  
#  # setting
#  n1 <- n2 <- 20
#  N <- c(n1,n2)
#  set.seed(123)
#  
#  # power
#  p.values <- matrix(NA,m,3)
#  for(i in 1:m){
#    x <- rnorm(n1,-1,2)
#    y <- rnorm(n2,0.5,1)
#    z <- c(x,y)
#    p.values[i,1] <- eqdist.nn(z,N,k)$p.value
#    p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
#    p.values[i,3] <- bd.test(x,y,num.permutations=R,seed=i*123)$p.value
#  }
#  alpha <- 0.05;
#  pow <- colMeans(p.values<alpha)
#  pow

## ----eval=FALSE---------------------------------------------------------------
#  k <- 3;R <- 500-1;m <- 200
#  
#  # setting
#  n1 <- n2 <- 20
#  N <- c(n1,n2)
#  set.seed(123)
#  
#  # power
#  p.values <- matrix(NA,m,3)
#  for(i in 1:m){
#    x <- rt(n1,df=1)
#    index<-sample(c(1,0),size=n2,replace=T,prob=c(0.5,0.5))
#    y<-index*rnorm(n1,0,1)+(1-index)*rnorm(n2,1,sqrt(0.5))
#    z <- c(x,y)
#    p.values[i,1] <- eqdist.nn(z,N,k)$p.value
#    p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
#    p.values[i,3] <- bd.test(x,y,num.permutations=R,seed=i*123)$p.value
#  }
#  alpha <- 0.05;
#  pow <- colMeans(p.values<alpha)
#  pow

## ----eval=FALSE---------------------------------------------------------------
#  k <- 3;R <- 500-1;m <- 200
#  
#  # setting
#  n1 <- 10
#  n2 <- 100
#  N <- c(n1,n2)
#  set.seed(123)
#  
#  # power
#  p.values <- matrix(NA,m,3)
#  for(i in 1:m){
#    x <- mvrnorm(n1,c(0,0),diag(1,2))
#    y <- mvrnorm(n2,c(0,1),diag(3,2))
#    z <- rbind(x,y)
#    p.values[i,1] <- eqdist.nn(z,N,k)$p.value
#    p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
#    p.values[i,3] <- bd.test(x,y,num.permutations=R,seed=i*123)$p.value
#  }
#  alpha <- 0.05;
#  pow <- colMeans(p.values<alpha)
#  pow

## -----------------------------------------------------------------------------
m <- 5000;set.seed(200);x <- numeric(m);u <- runif(m)
x[1] <- rnorm(1) # generate X_0
for(i in 2:m){
  xt <- x[i-1]
  y <- rnorm(1,mean = 0,sd = abs(xt))
  rt <- dcauchy(y)*dnorm(xt,0,abs(y))/(dcauchy(xt)*dnorm(y,0,abs(xt)))
  if(u[i]<=rt) x[i] <- y
  else x[i] <- xt 
}

## -----------------------------------------------------------------------------
data <- x[1001:m]
q <- seq(.1,.9,.1)
sample_deciles <- quantile(data,q)
cauchy_deciles <- qcauchy(q)

result <- rbind(sample_deciles,cauchy_deciles)
colnames(result) <- q
knitr::kable(result)

## -----------------------------------------------------------------------------
m <- 5000;set.seed(1234);x <- numeric(m);u <- runif(m)
x[1] <- rnorm(1) 
for(i in 2:m){
  xt <- x[i-1]
  y <- rnorm(1,mean = xt,sd = 1)
  rt <- dcauchy(y)*dnorm(xt,y,1)/(dcauchy(xt)*dnorm(y,xt,1))
  if(u[i]<=rt) x[i] <- y
  else x[i] <- xt 
}

data <- x[1001:m]
q <- seq(.1,.9,.1)
sample_deciles <- quantile(data,q)
cauchy_deciles <- qcauchy(q)

result <- rbind(sample_deciles,cauchy_deciles)
colnames(result) <- q
knitr::kable(result)

## -----------------------------------------------------------------------------
Gibbs_sample <- function(N,a,b,n){ # N is the length of chain
  X <- matrix(0,nrow = N,ncol = 2) 
  X[1,] <- c(floor(n/2),0.5)       # Initial value
  for(i in 2:N){
    y <- X[i-1,2]
    X[i,1] <- rbinom(1,size = n,prob = y) 
    x <- X[i,1]
    X[i,2] <- rbeta(1,x+a,n-x+b)          
  }
  return(X)
}

data <- Gibbs_sample(N=5000,a=5,b=5,n=50)
burn <- 1000
data <- data[(burn+1):5000,]        
plot(data,col="lightblue",pch=20,cex=0.6,xlab = "x",ylab = "y",ylim = c(0,1),xlim = c(0,50))

## ----eval=FALSE---------------------------------------------------------------
#  # Gelman-Rubin method
#  Gelman.Rubin <- function(psi){
#    psi <- as.matrix(psi)
#    n <- ncol(psi)
#    k <- nrow(psi)
#    psi.means <- rowMeans(psi)
#    B <- n*var(psi.means)
#    psi.w <- apply(psi,1,"var") # within variances
#    W <- mean(psi.w)
#    v.hat <- W*(n-1)/n + (B/n)
#    r.hat <- v.hat/W   # G-R statistic
#    return(r.hat)
#  }

## ----eval=FALSE---------------------------------------------------------------
#  # 9.3
#  cauchy.chain <- function(N,X1){
#    x <- rep(0,N)
#    x[1] <- X1
#    u <- runif(N)
#  for(i in 2:N){
#    xt <- x[i-1]
#    y <- rnorm(1,mean = xt,sd = 1)
#    rt <- dcauchy(y)*dnorm(xt,y,1)/(dcauchy(xt)*dnorm(y,xt,1))
#    if(u[i]<=rt) x[i] <- y
#    else x[i] <- xt
#  }
#    return(x)
#  }
#  
#  k <- 4
#  n <- 20000
#  b <- 2000  # burn-in length
#  
#  init <- c(-1.2,-1,-1.1,-0.5);set.seed(200) # Initial value
#  X <- matrix(0,nrow = k,ncol = n)
#  for(i in 1:k){
#    X[i,] <- cauchy.chain(n,init[i])
#  }
#  
#  # Calculate diagnostic statistics
#  psi <- t(apply(X,1,cumsum))
#  for(i in 1:nrow(psi)){
#    psi[i,] <- psi[i,]/(1:ncol(psi))
#  }
#  print(Gelman.Rubin(psi))
#  
#  # plot
#  rhat <- rep(0:n)
#  for(j in (b+1):n){
#    rhat[j] <- Gelman.Rubin(psi[,1:j])
#  }
#  plot(rhat[(b+1):n],type="l",xlab="",ylab="R")
#  abline(h=1.2,lty=2)

## ----eval=FALSE---------------------------------------------------------------
#  # 9.8
#  k <- 4
#  n <- 15000
#  burn <- 1000  # burn-in length
#  
#  
#  X <- matrix(0,nrow = k,ncol = n)
#  Y <- matrix(0,nrow = k,ncol = n)
#  set.seed(12345)
#  for(i in 1:k){
#    chain <- Gibbs_sample(n,5,5,30)
#    X[i,] <- chain[,1]
#    Y[i,] <- chain[,2]
#  }
#  
#  # Calculate diagnostic statistics of X
#  psi <- t(apply(X,1,cumsum))
#  for(i in 1:nrow(psi)){
#    psi[i,] <- psi[i,]/(1:ncol(psi))
#  }
#  print(Gelman.Rubin(psi))
#  
#  rhat <- rep(0:n)
#  for(j in (b+1):n){
#    rhat[j] <- Gelman.Rubin(psi[,1:j])
#  }
#  plot(rhat[(b+1):n],type="l",xlab="X",ylab="R")
#  abline(h=1.2,lty=2)
#  
#  # Calculate diagnostic statistics of Y
#  psi <- t(apply(Y,1,cumsum))
#  for(i in 1:nrow(psi)){
#    psi[i,] <- psi[i,]/(1:ncol(psi))
#  }
#  print(Gelman.Rubin(psi))
#  
#  rhat <- rep(0:n)
#  for(j in (b+1):n){
#    rhat[j] <- Gelman.Rubin(psi[,1:j])
#  }
#  plot(rhat[(b+1):n],type="l",xlab="Y",ylab="R")
#  abline(h=1.2,lty=2)

## -----------------------------------------------------------------------------
ak <- function(k,a){
  d <- length(a)      
  mode.a2 <- sum(a^2) 
  g1 <- exp(lgamma((d+1)/2)+lgamma(k+3/2)-lgamma(k+d/2+1)-lgamma(k+1)-k*log(2)) # k!=gamma(k+1)
  #g2 <- exp(-sum(log(1:k))-k*log(2))  # 1/(k!*2^k)
  return((-1)^k*mode.a2^(k+1)*g1/((2*k+1)*(2*k+2)))
}

## -----------------------------------------------------------------------------
summ <- function(K,a){
  sum(ak((0:K),a))
}

## -----------------------------------------------------------------------------
summ(200,c(1,2))

## -----------------------------------------------------------------------------
# 11.4
k <- c(4:25,100,500,1000)

s <- function(k,a){
  # S
  return(pt(sqrt(a^2*k/(k+1-a^2)),df=k,lower.tail = F))
}

delta.s <- function(k,a){
  s(k,a)-s(k-1,a)
}

n <- length(k)
A <- numeric(n)
for(i in 1:n){
  A[i] <- uniroot(delta.s,interval=c(1,2),k=k[i])$root
}

result <- cbind(k,A)
knitr::kable(result)

## -----------------------------------------------------------------------------
# 11.5
f1 <- function(k,a){
  integrate(dt, lower = 0, upper = sqrt(a^2*k/(k+1-a^2)) ,df = k)$value

}
f <- function(k, a){
  f1(k, a)-f1(k-1, a)
}

k <- c(4:25, 100, 500, 1000)
n <- length(k)
a <- numeric(n)
for(i in 1:n){
 a[i] <- uniroot(f, interval=c(1, 2), k = k[i])$root
}

result <- cbind(k,a)
knitr::kable(result)

## -----------------------------------------------------------------------------
y <- c(0.54,0.48,0.33,0.43,1.00,1.00,0.91,1.00,0.21,0.85)

f <- function(lambda){
  10/(6.75+3/lambda)
}

lambda.EM <- function(N,epsilon){
  lambda0 <- 1 
  for(i in 1:N){
    lambda1 <- f(lambda0)
    if(abs(lambda1-lambda0) < epsilon){
      return(lambda1) 
    }
    else lambda0 <- lambda1
  }
  print("EM algorithm is not covergent")
}

N <- 1000
epsilon <- 10^(-8)
EM <- lambda.EM(N,epsilon)

## -----------------------------------------------------------------------------
lambda.MLE <- 7/6.75
print(c(EM=EM,MLE=lambda.MLE))

## ----eval=FALSE---------------------------------------------------------------
#  trims <- c(0,0.1,0.2,0.5)
#  x <- rcauchy(100)
#  
#  lapply(trims,function(trim) mean(x,trim=trim))
#  lapply(trims, mean,x=x)

## -----------------------------------------------------------------------------
rsq <- function(mod) summary(mod)$r.squared

# 3
formulas <- list(
  mpg ~ disp,
  mpg ~ I(1/disp),
  mpg ~ disp + wt,
  mpg ~ I(1/disp) + wt
)
model1 <- lapply(formulas,lm,data=mtcars) # Use lapply()
unlist(lapply(model1,rsq)) # R2

loop1 <- vector("list",length = length(formulas)) # Use loops
for(i in seq_along(formulas)){
  loop1[[i]] <- lm(formulas[[i]],data = mtcars)
}
unlist(lapply(loop1,rsq)) # R2

# 4
bootstaps <- lapply(1:10, function(i){
  rows <- sample(1:nrow(mtcars),replace = T) # Resampling
  mtcars[rows,]
})

model2 <- lapply(bootstaps,lm,formula=mpg ~ disp) # Use lapply()
unlist(lapply(model2,rsq)) # R2

loop2 <- vector("list",length = length(bootstaps))
for(i in seq_along(bootstaps)){
  loop2[[i]] <- lm(mpg ~ disp,data=bootstaps[[i]])
}
unlist(lapply(loop2,rsq)) # R2

## -----------------------------------------------------------------------------
# numeric data frame
vapply(mtcars, sd ,numeric(1))  # Specify the return value type as numeric

## -----------------------------------------------------------------------------
# mixed data frame
data("iris")
sapply(iris, class)  
vapply(iris[vapply(iris, is.numeric, logical(1))],sd,numeric(1))

## ----eval=FALSE---------------------------------------------------------------
#  library(parallel)
#  detectCores()  # Check the number of cores currently available
#  
#  mcsapply <- function(x,fun){
#    cl <- makeCluster(detectCores())
#    results <- parSapply(cl,x,fun)
#    stopCluster(cl)
#    return(results)
#  }
#  
#  boot_lm <- function(i){
#    boot_df <- function(x) x[sample(nrow(x),replace = T),]
#    summary(lm(mpg~disp,data=boot_df(mtcars)))$r.squared
#  }
#  
#  n <- 4000
#  system.time(sapply(1:n,boot_lm))
#  system.time(mcsapply(1:n,boot_lm))

## ----eval=FALSE---------------------------------------------------------------
#  #include <Rcpp.h>
#  using namespace Rcpp;
#  
#  // [[Rcpp::export]]
#  NumericMatrix gibbsC(int N, double a, double b, int n, int thin) {
#    NumericMatrix mat(N, 2);
#    double x = floor(n/2), y = 0.5;
#    for(int i = 0; i < N; i++) {
#      for(int j = 0; j < thin; j++) {
#        x = rbinom(1, n, y )[0];
#        y = rbeta(1, x+a, n-x+b)[0];
#      }
#      mat(i, 0) = x;
#      mat(i, 1) = y;
#    }
#    return(mat);
#  }

## ----eval=FALSE---------------------------------------------------------------
#  library(Rcpp)
#  dir_cpp <- '../Rcpp/'
#  sourceCpp(paste0(dir_cpp,"gibbsC.cpp"))

## ----eval=FALSE---------------------------------------------------------------
#  gibbsR <- function(N,a,b,n, thin){ # N is the length of chain
#    X <- matrix(0,nrow = N,ncol = 2)
#    x <- floor(n/2)
#    y <- 0.5
#    for(i in 2:N){
#      for(j in 1:thin){
#        x <- rbinom(1,size = n,prob = y) #
#        y <- rbeta(1,x+a,n-x+b)
#      }
#      X[i,] <- c(x,y)
#    }
#    return(X)
#  }

