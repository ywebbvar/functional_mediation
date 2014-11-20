---
title: "Simulation scalar-scalar-function"
author: "Yenny Webb-Vargas"
date: "Wednesday, September 17, 2014"
output: html_document
---

This creates a simulation for a scalar-scalar-function mediation set-up. It is based on the "Simulation using splines. Original. In Gypsie", in file `Simulation_model_based_fRegress.Rmd`.


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
if(cluster){
  source('/home/bst/student/ywebbvar/Github/functional_mediation/ssf_Mediation.R')
  }else{
  source('~/GitHub/functional_mediation/ssf_Mediation.R')  
  }
```

Defining bootstrap function:

```r
fbootstrap_ML <- function(dta, index){
  dta = dta[index,]
  M = dta[,"M"]
  X = dta[,"X"]
  Y = t(dta[,-c(grep("M",colnames(dta)), grep("X",colnames(dta)))])
  
  result = ssf_Mediation(X,Y,M,outcomeMethod="fosr2s", nbasis=nbasis,norder=norder,plot=FALSE, boot=TRUE)
  return(result)
}
```


```r
sub=20
nbasis = 50
norder = 6
lambda = 1e-8

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
sd_eta     = 1
```

Create parameters: covariate and treatment effects

```r
# Create delta2
set.seed(73390)
delta2_coef = c(rep(0,5), runif(5,2,5), runif(12,-2,3),rep(0,8)) 
delta2      = eval.bspline%*%delta2_coef

# Create gamma
set.seed(2939)
gamma_coef = c(rep(0,11), runif(3,3,5),runif(3,6,8),runif(3,1,4), rep(0,10)) 
gamma      = eval.bspline%*%gamma_coef

# Create beta
set.seed(38923)
beta_coef = c(rep(0,15), runif(5,3,10), rep(0,10)) 
beta      = eval.bspline%*%beta_coef

delta1 = 5
alpha  = 3

if(plots) matplot(cbind(delta2,gamma,beta), type="l")
```

Setting seed:

```r
seed = 279500
set.seed(seed)
```


### Simulation 4
In this simulation, $M$ depends on $X$, $Y$ depends on $M$ and $X$.


```r
Pall  = matrix(0, iters, len)
Pall2 = matrix(0, iters, len)
d2    = matrix(0, iters, len)
B     = matrix(0, iters, len)
G     = matrix(0, iters, len)
d1    = rep(0, iters)
A     = rep(0, iters)

for(i in 1:iters){
  X          = rbinom(sub,1,0.5)
  
  epsilon    = rnorm(sub, 0, sd_epsilon) 

  design_M = cbind(int=1,X=X)
  coeff_M  = cbind(delta1,alpha)
  colnames(coeff_M) = c("delta1","alpha")

  M    = t(tcrossprod(coeff_M, design_M))+epsilon
  M    = as.numeric(M)
 
  eta  = matrix(rnorm(sub*len,0,sd_eta), ncol=len,nrow=(sub))

  design_M_Y = cbind(int=1,X=X,M=M)
  coeff_M_Y  = cbind(delta2,gamma, beta)
  colnames(coeff_M_Y) = c("delta2","gamma","beta")


  Y    <- t(tcrossprod(coeff_M_Y, design_M_Y)) + eta
  

  dta  = cbind('Y'=Y,'X'=X,'M' = M)
  
  out  = ssf_Mediation(X,t(Y),M,outcomeMethod="fosr2s", nbasis=nbasis,norder=norder,plot=plots, boot=FALSE)
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
d2[i,] = out$y_path$delta2_function
B[i,]  = out$y_path$beta_function
G[i,]  = out$y_path$gamma_function
d1[i]  = out$m_path$delta1
A[i]   = out$m_path$alpha
  
  if(plots)matplot(t(stat$t), type="l", col = rgb(200,1,1, 30, maxColorValue = 255))
  }
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9, 3, 1 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5, 7, 9, 8, 4, 1, 10 encountered errors in user code, all
## values of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 8, 7 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 3, 8, 7, 6, 10 encountered errors in user code, all values
## of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 8, 4, 3, 2 encountered errors in user code, all values of
## the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 8, 6 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 6, 7, 9 encountered errors in user code, all values of the
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
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 3, 9 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2, 1, 8, 7, 3 encountered errors in user code, all values
## of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9, 2, 7 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus): all
## scheduled cores encountered errors in user code
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 10, 6, 4, 8, 2 encountered errors in user code, all values
## of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5, 9, 3, 1 encountered errors in user code, all values of
## the jobs will be affected
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
## scheduled cores 10 encountered errors in user code, all values of the jobs
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
## scheduled cores 9, 2, 4, 6 encountered errors in user code, all values of
## the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 7 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 6, 10, 1 encountered errors in user code, all values of
## the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2 encountered errors in user code, all values of the jobs
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
## scheduled cores 7 encountered errors in user code, all values of the jobs
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
## scheduled cores 5, 4, 7 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 6, 8 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1 encountered error in user code, all values of the job
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4 encountered errors in user code, all values of the jobs
## will be affected
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
## scheduled cores 6 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 3, 7 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9, 2, 8, 1, 5, 4, 7, 3 encountered errors in user code,
## all values of the jobs will be affected
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
## scheduled cores 6 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 8 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2, 7 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 7 encountered errors in user code, all values of the jobs
## will be affected
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
## scheduled cores 10 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 3, 1 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1, 7, 3 encountered error in user code, all values of the
## job will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 6 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 6, 3 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 3 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 10, 2, 9, 5, 4 encountered errors in user code, all values
## of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 6, 1, 4, 7, 10, 3, 9, 8 encountered errors in user code,
## all values of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 6, 5, 4, 7 encountered errors in user code, all values of
## the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 3, 10 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus): all
## scheduled cores encountered errors in user code
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4, 7, 10 encountered errors in user code, all values of
## the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 3, 5 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus): all
## scheduled cores encountered errors in user code
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9 encountered errors in user code, all values of the jobs
## will be affected
```

```r
save(d2,B,G,d1,A,Pall,Pall2, file=paste0("SSF_Simulation_4_",format(Sys.time(), "%Y%m%d-%H%M"),".Rdata"))
```


### Null Simulation 1. 
In this simulation, $M$ is pure noise (and is not dependent on $X$), and $Y$ only depends on $X$.


```r
Pall  = matrix(0, iters, len)
Pall2 = matrix(0, iters, len)
d2    = matrix(0, iters, len)
B     = matrix(0, iters, len)
G     = matrix(0, iters, len)
d1    = rep(0, iters)
A     = rep(0, iters)

for(i in 1:iters){
  X          = rbinom(sub,1,0.5)
  
  epsilon    = rnorm(sub, 0, sd_epsilon) 

  design_M = cbind(int=1,X=X)
  coeff_M  = cbind(0,0)
  colnames(coeff_M) = c("delta1","alpha")

  M    = t(tcrossprod(coeff_M, design_M))+epsilon
  M    = as.numeric(M)
 
  eta  = matrix(rnorm(sub*len,0,sd_eta), ncol=len,nrow=(sub))

  design_M_Y = cbind(int=1,X=X,M=M)
  coeff_M_Y  = cbind(delta2,gamma, rep(0, length(beta)))
  colnames(coeff_M_Y) = c("delta2","gamma","beta")


  Y    <- t(tcrossprod(coeff_M_Y, design_M_Y)) + eta

  dta  = cbind('Y'=Y,'X'=X,'M' = M)
  
  out  = ssf_Mediation(X,t(Y),M,outcomeMethod="fosr2s", nbasis=nbasis,norder=norder,plot=plots, boot=FALSE)
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
d2[i,] = out$y_path$delta2_function
B[i,]  = out$y_path$beta_function
G[i,]  = out$y_path$gamma_function
d1[i]  = out$m_path$delta1
A[i]   = out$m_path$alpha
  
  if(plots)matplot(t(stat$t), type="l", col = rgb(200,1,1, 30, maxColorValue = 255))
  }
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1 encountered error in user code, all values of the job
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2, 7 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5, 6, 8 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 6 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2, 1 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1, 9, 7, 10, 6, 4, 8, 3, 2 encountered error in user code,
## all values of the job will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2, 6 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 6, 2 encountered errors in user code, all values of the
## jobs will be affected
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
## scheduled cores 5 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2, 10, 9 encountered errors in user code, all values of
## the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5, 2 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 6 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 7, 4 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 6, 4 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9 encountered errors in user code, all values of the jobs
## will be affected
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
## scheduled cores 7, 8 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1, 6, 8, 5 encountered error in user code, all values of
## the job will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 10 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9, 8 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4, 10, 8, 2, 9 encountered errors in user code, all values
## of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5, 1, 8 encountered errors in user code, all values of the
## jobs will be affected
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
## scheduled cores 9, 10, 4, 1, 6, 2, 3 encountered errors in user code, all
## values of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus): all
## scheduled cores encountered errors in user code
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1, 8 encountered error in user code, all values of the job
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 6 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 8, 2 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4, 2 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 8, 9 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 7 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 10, 2, 8, 3, 6, 4, 5 encountered errors in user code, all
## values of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 8, 3, 2 encountered errors in user code, all values of the
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
## scheduled core 1 encountered error in user code, all values of the job
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1, 9, 2, 4, 5, 3, 7 encountered error in user code, all
## values of the job will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 7, 3, 2, 5, 10, 1, 9, 4 encountered errors in user code,
## all values of the jobs will be affected
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
## scheduled cores 7, 9, 2 encountered errors in user code, all values of the
## jobs will be affected
```

```r
save(d2,B,G,d1,A,Pall,Pall2, file=paste0("SSF_Simulation_1_",format(Sys.time(), "%Y%m%d-%H%M"),".Rdata"))
```



### Null Simulation 2. 
In this simulation, $M$ depends on $X$, but $Y$ only depends on $X$ and not on $M$.


```r
Pall  = matrix(0, iters, len)
Pall2 = matrix(0, iters, len)
d2    = matrix(0, iters, len)
B     = matrix(0, iters, len)
G     = matrix(0, iters, len)
d1    = rep(0, iters)
A     = rep(0, iters)

for(i in 1:iters){
  X          = rbinom(sub,1,0.5)
  
  epsilon    = rnorm(sub, 0, sd_epsilon) 

  design_M = cbind(int=1,X=X)
  coeff_M  = cbind(delta1,alpha)
  colnames(coeff_M) = c("delta1","alpha")

  M    = t(tcrossprod(coeff_M, design_M))+epsilon
  M    = as.numeric(M)
 
  eta  = matrix(rnorm(sub*len,0,sd_eta), ncol=len,nrow=(sub))

  design_M_Y = cbind(int=1,X=X,M=M)
  coeff_M_Y  = cbind(delta2,gamma, rep(0, length(beta)))
  colnames(coeff_M_Y) = c("delta2","gamma","beta")


  Y    <- t(tcrossprod(coeff_M_Y, design_M_Y)) + eta

  dta  = cbind('Y'=Y,'X'=X,'M' = M)
  
  out  = ssf_Mediation(X,t(Y),M,outcomeMethod="fosr2s", nbasis=nbasis,norder=norder,plot=plots, boot=FALSE)
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
d2[i,] = out$y_path$delta2_function
B[i,]  = out$y_path$beta_function
G[i,]  = out$y_path$gamma_function
d1[i]  = out$m_path$delta1
A[i]   = out$m_path$alpha
  
  if(plots)matplot(t(stat$t), type="l", col = rgb(200,1,1, 30, maxColorValue = 255))
  }
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 8, 9, 1, 3, 7, 2, 6 encountered errors in user code, all
## values of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1, 4 encountered error in user code, all values of the job
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 8, 1, 2 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 3 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 10, 7, 2 encountered errors in user code, all values of
## the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 10 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1, 2 encountered error in user code, all values of the job
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 7, 2 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 6, 2 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 10, 9, 2, 3 encountered errors in user code, all values of
## the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 10, 3, 1 encountered errors in user code, all values of
## the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 8 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 3 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4, 9, 6, 7, 2 encountered errors in user code, all values
## of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5, 9, 10, 4, 2, 6 encountered errors in user code, all
## values of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 3, 2 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9, 10, 8, 7 encountered errors in user code, all values of
## the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9, 2, 1, 7 encountered errors in user code, all values of
## the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 7, 9, 3 encountered errors in user code, all values of the
## jobs will be affected
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
## scheduled cores 7 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1, 9, 10, 6, 5, 4, 2 encountered error in user code, all
## values of the job will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5, 2 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 7, 6, 4, 9, 5, 3 encountered errors in user code, all
## values of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 8, 7, 2, 6 encountered errors in user code, all values of
## the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 10 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 7 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 8, 3, 10 encountered errors in user code, all values of
## the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 3 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 8 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9, 5 encountered errors in user code, all values of the
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
## scheduled cores 9, 6, 2, 1, 5, 7, 10, 8 encountered errors in user code,
## all values of the jobs will be affected
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
## scheduled cores 9, 7 encountered errors in user code, all values of the
## jobs will be affected
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
## scheduled cores 7 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 8, 1, 4, 5, 7 encountered errors in user code, all values
## of the jobs will be affected
```

```r
save(d2,B,G,d1,A,Pall,Pall2, file=paste0("SSF_Simulation_2_",format(Sys.time(), "%Y%m%d-%H%M"),".Rdata"))
```


### Null simulation 3
In this simulation, $M$ depends on $\delta_1$ for all subjects (and does not depend on $X$), $Y$ only depends on $M$ and not on $X$.


```r
Pall  = matrix(0, iters, len)
Pall2 = matrix(0, iters, len)
d2    = matrix(0, iters, len)
B     = matrix(0, iters, len)
G     = matrix(0, iters, len)
d1    = rep(0, iters)
A     = rep(0, iters)

for(i in 1:iters){
  X          = rbinom(sub,1,0.5)
  
  epsilon    = rnorm(sub, 0, sd_epsilon) 

  design_M = cbind(int=1,X=X)
  coeff_M  = cbind(delta1,0)
  colnames(coeff_M) = c("delta1","alpha")

  M    = t(tcrossprod(coeff_M, design_M))+epsilon
  M    = as.numeric(M)
 
  eta  = matrix(rnorm(sub*len,0,sd_eta), ncol=len,nrow=(sub))

  design_M_Y = cbind(int=1,X=X,M=M)
  coeff_M_Y  = cbind(delta2,rep(0, length(gamma)), beta)
  colnames(coeff_M_Y) = c("delta2","gamma","beta")


  Y    <- t(tcrossprod(coeff_M_Y, design_M_Y)) + eta

  dta  = cbind('Y'=Y,'X'=X,'M' = M)
  
  out  = ssf_Mediation(X,t(Y),M,outcomeMethod="fosr2s", nbasis=nbasis,norder=norder,plot=plots, boot=FALSE)
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
d2[i,] = out$y_path$delta2_function
B[i,]  = out$y_path$beta_function
G[i,]  = out$y_path$gamma_function
d1[i]  = out$m_path$delta1
A[i]   = out$m_path$alpha
  
  if(plots)matplot(t(stat$t), type="l", col = rgb(200,1,1, 30, maxColorValue = 255))
  }
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
## scheduled cores 9 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus): all
## scheduled cores encountered errors in user code
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1 encountered error in user code, all values of the job
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 8, 7 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 7 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 8 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4, 7 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 3, 2, 8, 6, 1 encountered errors in user code, all values
## of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 10 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 10, 3, 8, 4, 9, 2, 1, 6, 7 encountered errors in user
## code, all values of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1, 3 encountered error in user code, all values of the job
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 6 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2, 4, 3, 8, 7, 9 encountered errors in user code, all
## values of the jobs will be affected
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
## scheduled cores 2 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 7, 10, 3, 6 encountered errors in user code, all values of
## the jobs will be affected
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
## scheduled cores 10, 5 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9, 5, 10, 7 encountered errors in user code, all values of
## the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 7, 8, 9 encountered errors in user code, all values of the
## jobs will be affected
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

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 10 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 7, 6, 9, 4, 10, 1, 3, 2 encountered errors in user code,
## all values of the jobs will be affected
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
## scheduled cores 9, 2, 5 encountered errors in user code, all values of the
## jobs will be affected
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
## scheduled cores 9 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 3, 2 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9, 2 encountered errors in user code, all values of the
## jobs will be affected
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
## scheduled cores 7 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 8 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 3, 7, 4 encountered errors in user code, all values of the
## jobs will be affected
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
## scheduled cores 8 encountered errors in user code, all values of the jobs
## will be affected
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
## scheduled cores 5 encountered errors in user code, all values of the jobs
## will be affected
```

```r
save(d2,B,G,d1,A,Pall,Pall2, file=paste0("SSF_Simulation_3_",format(Sys.time(), "%Y%m%d-%H%M"),".Rdata"))
```
