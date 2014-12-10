---
title: "Simulation scalar-function-scalar"
author: "Yenny Webb-Vargas"
date: "Monday, December 8, 2014"
output: html_document
---

This implements simulation for testing my functional mediation function. It uses new function for mediation models scalar-function-scalar that implements gam estimation

Load the library 'refund'. It also loads the library 'fda' among many others.


```r
cluster = TRUE
plots   = FALSE
```


```r
library(refund)
```

```
## Loading required package: fda
## Loading required package: splines
## Loading required package: Matrix
## Loading required package: methods
## 
## Attaching package: 'fda'
## 
## The following object is masked from 'package:graphics':
## 
##     matplot
```

```r
library(boot)
```

```
## 
## Attaching package: 'boot'
## 
## The following object is masked from 'package:fda':
## 
##     melanoma
```

```r
library(mgcv)
```

```
## Loading required package: nlme
## This is mgcv 1.8-3. For overview type 'help("mgcv-package")'.
```

```r
if(cluster){
  source('/home/bst/student/ywebbvar/Github/functional_mediation/sfs_Mediation.R')
  source('/home/bst/student/ywebbvar/Github/functional_mediation/est_se_fgam.R')
  }else{
  source('~/GitHub/functional_mediation/sfs_Mediation.R')  
  source('~/GitHub/functional_mediation/est_se_fgam.R')  
  }
```

Defining bootstrap function:

```r
fbootstrap_ML <- function(dta, index){
  dta = dta[index,]
  Y = dta[,"Y"]
  X = dta[,"X"]
  M = t(dta[,seq(2+1,ncol(dta))])
  
  result = sfs_Mediation(X,Y,M,mediatorMethod="fosr2s", outcomeMethod="fgam",nbasis=nbasis,norder=norder,lambda=lambda, plot=FALSE, boot=TRUE)
  return(result)
}
```


```r
sub=20
nbasis = 50
norder = 6
lambda = 1e-10 #1e-14 when 0.01 sd_epsilon

if(cluster){
  iters   = 500  
  nboot   = 1000 
}else{
  iters   = 1  
  nboot   = 50 
}

len     = 50
T_sup   = 1
timevec = seq(0,T_sup, length.out=len)
  
# Create bspline basis set
basis        = create.bspline.basis(rangeval = c(0,T_sup), nbasis = 30, norder=6)
eval.bspline = eval.basis(seq(0,1,length.out=len), basis)

# Parameters for random errors
sd_epsilon = 1 
sd_eta     = 0.1 
```

Create parameters: covariate and treatment effects

```r
# Create delta
set.seed(2903)
delta1_coef = c(rep(0,5), runif(5,-2,0), runif(10,0,2),rep(0,10)) 
delta1      = eval.bspline%*%delta1_coef

# Create alpha
set.seed(493)
alpha_coef = c(rep(0,10), runif(10,1,4), rep(0,10)) 
alpha      = eval.bspline%*%alpha_coef

# Create beta
set.seed(84930)
beta_coef = c(rep(0,15), runif(5,3,10), rep(0,10)) 
beta      = eval.bspline%*%beta_coef

delta2 = 5
gamma  = 2

if(plots) matplot(cbind(delta1,alpha,beta), type="l")
```

Setting seed:

```r
seed = 632020
set.seed(seed)
```


### Simulation 4
In this simulation, $M$ depends on $Z$ and $X$, $Y$ only depends on $M$ and $X$.


```r
Pall  = matrix(0, iters, len)
Pall2 = matrix(0, iters, len)
Flag  = rep(0, iters)
A     = matrix(0, iters, len)
B     = matrix(0, iters, len)
CP    = rep(0, iters)
d2    = rep(0, iters)
ptime = list()

for(i in 1:iters){
  X          = rbinom(sub,1,0.5)
  
  epsilon    = matrix(rnorm(sub*len,0,sd_epsilon), ncol=len,nrow=(sub))

  design_M = cbind(int=1,X=X)
  coeff_M  = cbind(delta1,alpha)

  colnames(coeff_M) = c("delta1","alpha")

  M   = t(tcrossprod(coeff_M, design_M))+epsilon

  eta = rnorm(sub, 0, sd_eta) 

  design_M_Y = cbind(int=1,X=X,M=M)
  coeff_M_Y  = cbind(delta2,gamma, t(beta)*1/(ncol(M)-1))

  Y    = t(tcrossprod(coeff_M_Y, design_M_Y)) + eta
  Y    = as.numeric(Y)

  dta  = cbind('Y'=Y,'X'=X,'M' = M)
  
  out  = sfs_Mediation(X,Y,t(M),mediatorMethod="fosr2s", outcomeMethod="fgam",nbasis=nbasis,norder=norder,lambda=lambda, plot=plots, boot=FALSE)
  if(cluster==TRUE){
    stat = boot(data = dta, statistic = fbootstrap_ML,R = nboot,sim = "ordinary",stype = "i",parallel = "multicore",ncpus=10)
    }else{
    stat = boot(data = dta, statistic = fbootstrap_ML,R = nboot,sim = "ordinary",stype = "i")
    }

                
  
  # Calculating p-value
  Pall[i,] = 2*apply(cbind(colMeans(stat$t < 2*matrix(stat$t0, ncol=ncol(stat$t), nrow=nrow(stat$t), byrow = TRUE)),
                           colMeans(stat$t > 2*matrix(stat$t0, ncol=ncol(stat$t), nrow=nrow(stat$t), byrow = TRUE))),1, min)
  Pall2[i,] = 2*apply(cbind(colMeans(stat$t < matrix(stat$t0, ncol=ncol(stat$t), nrow=nrow(stat$t), byrow = TRUE)),
                            colMeans(stat$t > matrix(stat$t0, ncol=ncol(stat$t), nrow=nrow(stat$t), byrow = TRUE))), 1, min)
  A[i,] = out$afunction
  B[i,] = out$bfunction
  CP[i] = out$cp
  d2[i] = out$Y_intercept
  
  if(plots)matplot(t(stat$t), type="l", col = rgb(200,1,1, 30, maxColorValue = 255))
  }
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 3 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 7 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5, 3, 8, 10, 6, 2, 9 encountered errors in user code, all
## values of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1 encountered error in user code, all values of the job
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9, 10, 7 encountered errors in user code, all values of
## the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 8 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 6, 4, 1 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 6 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2, 9 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 8 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1 encountered error in user code, all values of the job
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 10 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4, 1 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 10, 9, 5, 2 encountered errors in user code, all values of
## the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1 encountered error in user code, all values of the job
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 10, 5, 1, 9, 3, 8, 4, 2, 6 encountered errors in user
## code, all values of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 7 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 3 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4, 7 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2, 1 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 6, 3 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2, 3 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 6 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 10 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 6, 2, 5, 7, 8, 9, 4, 1 encountered errors in user code,
## all values of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1, 5 encountered error in user code, all values of the job
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9, 1 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 8, 3, 6 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2, 6, 3 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5, 6, 4 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9, 7, 8, 1, 10, 2, 5, 4 encountered errors in user code,
## all values of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4, 9 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 3, 8 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5, 10 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 6 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1 encountered error in user code, all values of the job
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 10, 4 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 7 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 6 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 3 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9, 8 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5, 6, 10, 9, 8, 3 encountered errors in user code, all
## values of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9, 2, 7, 6, 3, 4 encountered errors in user code, all
## values of the jobs will be affected
```

```r
save(A,B,d2,CP,Pall,Pall2, file=paste0("Simulation_4_lesspenal_",format(Sys.time(), "%Y%m%d-%H%M"),".Rdata"))
```


### Null Simulation 1. 
In this simulation, $M$ is pure noise (and is not dependent on $X$), and $Y$ only depends on $X$.


```r
Pall  = matrix(0, iters, len)
Pall2 = matrix(0, iters, len)
Flag  = rep(0, iters)
A     = matrix(0, iters, len)
B     = matrix(0, iters, len)
CP    = rep(0, iters)
d2    = rep(0, iters)
ptime = list()

for(i in 1:iters){
  X          = rbinom(sub,1,0.5)
  
  epsilon    = matrix(rnorm(sub*len,0,sd_epsilon), ncol=len,nrow=(sub))

  M    = epsilon

  eta      <- rnorm(sub, 0, sd_eta) 

  design_M_Y <- cbind(int=1,X=X)
  coeff_M_Y  <- cbind(delta2,gamma)

  Y    <- t(tcrossprod(coeff_M_Y, design_M_Y)) + eta
  Y    <- as.numeric(Y)

  dta  = cbind('Y'=Y,'X'=X,'M' = M)
  
  out  = sfs_Mediation(X,Y,t(M),mediatorMethod="fosr2s", nbasis=nbasis,norder=norder,lambda=lambda, plot=plots, boot=FALSE)
  if(cluster==TRUE){
    stat = boot(data = dta, statistic = fbootstrap_ML,R = nboot,sim = "ordinary",stype = "i",parallel = "multicore",ncpus=10)
    }else{
    stat = boot(data = dta, statistic = fbootstrap_ML,R = nboot,sim = "ordinary",stype = "i")
    }

                
  
  # Calculating p-value
  Pall[i,] = 2*apply(cbind(colMeans(stat$t < 2*matrix(stat$t0, ncol=ncol(stat$t), nrow=nrow(stat$t), byrow = TRUE)),
                           colMeans(stat$t > 2*matrix(stat$t0, ncol=ncol(stat$t), nrow=nrow(stat$t), byrow = TRUE))),1, min)
  Pall2[i,] = 2*apply(cbind(colMeans(stat$t < matrix(stat$t0, ncol=ncol(stat$t), nrow=nrow(stat$t), byrow = TRUE)),
                            colMeans(stat$t > matrix(stat$t0, ncol=ncol(stat$t), nrow=nrow(stat$t), byrow = TRUE))), 1, min)
  A[i,] = out$afunction
  B[i,] = out$bfunction
  CP[i] = out$cp
  d2[i] = out$Y_intercept
  
  if(plots)matplot(t(stat$t), type="l", col = rgb(200,1,1, 30, maxColorValue = 255))
  }
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 6, 10, 4, 1, 8, 2 encountered errors in user code, all
## values of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 10 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1 encountered error in user code, all values of the job
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 8 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 8 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 8, 6, 10, 4, 7, 2 encountered errors in user code, all
## values of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5, 2 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 10 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1, 10 encountered error in user code, all values of the job
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1, 7 encountered error in user code, all values of the job
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2, 7, 4 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9, 5, 1, 2, 3, 6, 8, 4, 10 encountered errors in user
## code, all values of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 10, 5, 8, 7, 3, 2, 1 encountered errors in user code, all
## values of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1 encountered error in user code, all values of the job
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 10 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5, 2 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4, 10, 9, 2, 7, 3, 6, 8 encountered errors in user code,
## all values of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1 encountered error in user code, all values of the job
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 7, 5 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 8, 4, 9 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 6 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 7, 2 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 6 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1, 7, 5, 10 encountered error in user code, all values of
## the job will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9, 6, 2, 1, 4, 7, 8 encountered errors in user code, all
## values of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5, 1 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 7 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 10 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4, 6 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 3 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1 encountered error in user code, all values of the job
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 6 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 6 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 3, 6 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1, 2, 5, 10, 7, 6, 3, 9 encountered error in user code, all
## values of the job will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 3, 2 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 7, 4 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2, 4, 6, 5, 3 encountered errors in user code, all values
## of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 8 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 7 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5, 1, 6, 7 encountered errors in user code, all values of
## the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 8 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 7 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1 encountered error in user code, all values of the job
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 3, 9 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2, 7 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 7, 4, 10, 1, 9, 6, 3, 2 encountered errors in user code,
## all values of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 8 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 6 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 6, 3, 4, 10, 1, 9, 8, 5 encountered errors in user code,
## all values of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5, 7, 10, 3, 6, 1, 4 encountered errors in user code, all
## values of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 10 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 6 encountered errors in user code, all values of the jobs
## will be affected
```

```r
save(A,B,d2,CP,Pall,Pall2, file=paste0("Simulation_1_lesspenal_",format(Sys.time(), "%Y%m%d-%H%M"),".Rdata"))
```



### Null Simulation 2. 
In this simulation, $M$ depends on $X$, but $Y$ only depends on $X$ and not on $M$.


```r
Pall  = matrix(0, iters, len)
Pall2 = matrix(0, iters, len)
Flag  = rep(0, iters)
A     = matrix(0, iters, len)
B     = matrix(0, iters, len)
CP    = rep(0, iters)
d2    = rep(0, iters)
ptime = list()

for(i in 1:iters){
  X          = rbinom(sub,1,0.5)
  
  epsilon    = matrix(rnorm(sub*len,0,sd_epsilon), ncol=len,nrow=(sub))

  design_M = cbind(int=1,X=X)
  coeff_M  = cbind(delta1,alpha)
  colnames(coeff_M) = c("delta1","alpha")

  M    = t(tcrossprod(coeff_M, design_M))+epsilon

  eta      <- rnorm(sub, 0, sd_eta) 

  design_M_Y <- cbind(int=1,X=X)
  coeff_M_Y  <- cbind(delta2,gamma)

  Y    <- t(tcrossprod(coeff_M_Y, design_M_Y)) + eta
  Y    <- as.numeric(Y)

  dta  = cbind('Y'=Y,'X'=X,'M' = M)
  
  out  = sfs_Mediation(X,Y,t(M),mediatorMethod="fosr2s", nbasis=nbasis,norder=norder,lambda=lambda, plot=plots, boot=FALSE)
  if(cluster==TRUE){
    stat = boot(data = dta, statistic = fbootstrap_ML,R = nboot,sim = "ordinary",stype = "i",parallel = "multicore",ncpus=10)
    }else{
    stat = boot(data = dta, statistic = fbootstrap_ML,R = nboot,sim = "ordinary",stype = "i")
    }

                
  
  # Calculating p-value
  Pall[i,] = 2*apply(cbind(colMeans(stat$t < 2*matrix(stat$t0, ncol=ncol(stat$t), nrow=nrow(stat$t), byrow = TRUE)),
                           colMeans(stat$t > 2*matrix(stat$t0, ncol=ncol(stat$t), nrow=nrow(stat$t), byrow = TRUE))),1, min)
  Pall2[i,] = 2*apply(cbind(colMeans(stat$t < matrix(stat$t0, ncol=ncol(stat$t), nrow=nrow(stat$t), byrow = TRUE)),
                            colMeans(stat$t > matrix(stat$t0, ncol=ncol(stat$t), nrow=nrow(stat$t), byrow = TRUE))), 1, min)
  A[i,] = out$afunction
  B[i,] = out$bfunction
  CP[i] = out$cp
  d2[i] = out$Y_intercept
  
  if(plots)matplot(t(stat$t), type="l", col = rgb(200,1,1, 30, maxColorValue = 255))
  }
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 10, 8, 4, 3, 1, 7 encountered errors in user code, all
## values of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus): all
## scheduled cores encountered errors in user code
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2, 6 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 10, 8 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2, 3, 5, 1, 4 encountered errors in user code, all values
## of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 3, 6 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 3, 4, 7, 6 encountered errors in user code, all values of
## the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4, 8 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 6, 2 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 6 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 3 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 10, 1 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2, 3, 8, 9 encountered errors in user code, all values of
## the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 7 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5, 9 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 7 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 10 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus): all
## scheduled cores encountered errors in user code
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1 encountered error in user code, all values of the job
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2, 1 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 10 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 8 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 7, 8 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1 encountered error in user code, all values of the job
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus): all
## scheduled cores encountered errors in user code
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1 encountered error in user code, all values of the job
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 8 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 8, 7, 2 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 7, 8, 6 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 10 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 10 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 10, 6 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1, 4, 7, 2, 6, 8, 9 encountered error in user code, all
## values of the job will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4, 8 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5, 2, 9, 1 encountered errors in user code, all values of
## the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 10, 9, 2, 5, 8, 6, 4, 3, 1 encountered errors in user
## code, all values of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 3 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2, 7, 5 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9, 2 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 10 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 6 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5 encountered errors in user code, all values of the jobs
## will be affected
```

```r
save(A,B,d2,CP,Pall,Pall2, file=paste0("Simulation_2_lesspenal_",format(Sys.time(), "%Y%m%d-%H%M"),".Rdata"))
```


### Null simulation 3
In this simulation, $M$ depends on $\delta_1$ for all subjects (and does not depend on $X$), $Y$ only depends on $M$ and not on $X$.


```r
Pall  = matrix(0, iters, len)
Pall2 = matrix(0, iters, len)
Flag  = rep(0, iters)
A     = matrix(0, iters, len)
B     = matrix(0, iters, len)
CP    = rep(0, iters)
d2    = rep(0, iters)

for(i in 1:iters){
  X          = rbinom(sub,1,0.5)
  
  epsilon    = matrix(rnorm(sub*len,0,sd_epsilon), ncol=len,nrow=(sub))

  design_M = rep(1,sub)
  coeff_M  = cbind(delta1)
  colnames(coeff_M) = c("delta1")

  M    = t(tcrossprod(coeff_M, design_M))+epsilon

  eta      <- rnorm(sub, 0, sd_eta) 

  design_M_Y <- cbind(int=1,X=X,M=M)
  coeff_M_Y  <- cbind(delta2,0, t(beta)*1/(ncol(M)-1))

  Y    <- t(tcrossprod(coeff_M_Y, design_M_Y)) + eta
  Y    <- as.numeric(Y)

  dta  = cbind('Y'=Y,'X'=X,'M' = M)
  
  out  = sfs_Mediation(X,Y,t(M),mediatorMethod="fosr2s", nbasis=nbasis,norder=norder,lambda=lambda, plot=plots, boot=FALSE)
  if(cluster==TRUE){
    stat = boot(data = dta, statistic = fbootstrap_ML,R = nboot,sim = "ordinary",stype = "i",parallel = "multicore",ncpus=10)
    }else{
    stat = boot(data = dta, statistic = fbootstrap_ML,R = nboot,sim = "ordinary",stype = "i")
    }

                
  
  # Calculating p-value
  Pall[i,] = 2*apply(cbind(colMeans(stat$t < 2*matrix(stat$t0, ncol=ncol(stat$t), nrow=nrow(stat$t), byrow = TRUE)),
                           colMeans(stat$t > 2*matrix(stat$t0, ncol=ncol(stat$t), nrow=nrow(stat$t), byrow = TRUE))),1, min)
  Pall2[i,] = 2*apply(cbind(colMeans(stat$t < matrix(stat$t0, ncol=ncol(stat$t), nrow=nrow(stat$t), byrow = TRUE)),
                            colMeans(stat$t > matrix(stat$t0, ncol=ncol(stat$t), nrow=nrow(stat$t), byrow = TRUE))), 1, min)
  A[i,] = out$afunction
  B[i,] = out$bfunction
  CP[i] = out$cp
  d2[i] = out$Y_intercept
  
  if(plots)matplot(t(stat$t), type="l", col = rgb(200,1,1, 30, maxColorValue = 255))
  }
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 10, 4, 5, 6, 8, 9, 1, 7 encountered errors in user code,
## all values of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 6 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 7, 4, 9, 5 encountered errors in user code, all values of
## the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 10, 8 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 10 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 10, 9, 4, 5, 8 encountered errors in user code, all values
## of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 6, 3 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2, 1, 3 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1, 2 encountered error in user code, all values of the job
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 7 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 7, 1, 6, 9, 10, 4, 8, 5, 2 encountered errors in user
## code, all values of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2, 3, 1 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2, 4 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1 encountered error in user code, all values of the job
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1 encountered error in user code, all values of the job
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4, 8, 7, 3 encountered errors in user code, all values of
## the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2, 7 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1 encountered error in user code, all values of the job
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 7 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1 encountered error in user code, all values of the job
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 7 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 6, 7, 2 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 8, 1 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 6 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 10 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1, 6 encountered error in user code, all values of the job
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9, 6, 10 encountered errors in user code, all values of
## the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 3, 7 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus): all
## scheduled cores encountered errors in user code
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2, 10 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 8 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1 encountered error in user code, all values of the job
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 6, 4, 3 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 3 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus): all
## scheduled cores encountered errors in user code
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus): all
## scheduled cores encountered errors in user code
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9, 10 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 7 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 3 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 3 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 8, 10 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 8 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 7 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 3 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2, 8 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1 encountered error in user code, all values of the job
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 8 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5, 3, 7, 4, 1, 2 encountered errors in user code, all
## values of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 6 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 6 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 3, 4 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 7, 10, 6, 9, 5, 2 encountered errors in user code, all
## values of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1 encountered error in user code, all values of the job
## will be affected
```

```r
save(A,B,d2,CP,Pall,Pall2, file=paste0("Simulation_3_lesspenal_",format(Sys.time(), "%Y%m%d-%H%M"),".Rdata"))
```
