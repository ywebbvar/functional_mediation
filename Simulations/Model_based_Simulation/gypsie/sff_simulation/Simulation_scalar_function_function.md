---
title: "Simulation using splines. Original. In Gypsie"
author: "Yenny Webb-Vargas"
date: "Monday, November 17, 2014"
output: html_document
---

This implements simulation for testing my functional mediation function. 

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
  source('/home/bst/student/ywebbvar/Github/functional_mediation/sff_Mediation.R')
  source('/home/bst/student/ywebbvar/Github/functional_mediation/est_se_fgam.R')
  source('/home/bst/student/ywebbvar/Github/functional_mediation/lplot.R')
  source('/home/bst/student/ywebbvar/Github/functional_mediation/sff_vec2list.R')
  }else{
  source('~/GitHub/functional_mediation/sff_Mediation.R') 
  source('~/GitHub/functional_mediation/est_se_fgam.R')
  source('~/GitHub/functional_mediation/lplot.R')
  source('~/GitHub/functional_mediation/sff_vec2list.R')
  }
```

Defining bootstrap function:

```r
create_functional_y = function(x, m, etai, delta2=delta2, gamma=gamma, beta=beta){
  m_mat = matrix(m, ncol=length(m), nrow=length(m), byrow=T)
  y = delta2 + gamma*x + rowSums(beta*m_mat*1/(ncol(M)-1)) + etai
  y
}

fbootstrap_ML = function(dta, index){
  dta = dta[index,]
  Y = dta[,grep("Y", names(dta))]
  X = dta[,grep("X", names(dta))]
  M = t(dta[,grep("M", names(dta))])
  
  result <- tryCatch(
    {sff_Mediation(X,Y,M,mediatorMethod="fosr2s",nbasis=nbasis,norder=norder,lambda=lambda, plot=FALSE, boot=TRUE)
     },
    error=function(cond){
      message(cond)
      return(rep(NA,(2*nrow(M)^2+5*nrow(M))))
    }
    )
  return(result)
}
```

Defining 


```r
sub=20
nbasis = 50
norder = 6
lambda = 1e-10

if(cluster){
  iters   = 500  
  nboot   = 1000 
}else{
  iters   = 1  
  nboot   = 100 
}

len     = 40
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

# Create delta2
set.seed(73390)
delta2_coef = c(rep(0,5), runif(5,2,5), runif(12,-2,3),rep(0,8)) 
delta2      = eval.bspline%*%delta2_coef

# Create gamma
set.seed(2939)
gamma_coef = c(rep(0,11), runif(3,3,5),runif(3,4,6),runif(3,1,4), rep(0,10)) 
gamma      = eval.bspline%*%gamma_coef

# Create beta
set.seed(84930)
beta_coef = c(rep(0,15), runif(10,6,8), runif(5,2,3)) 
beta_line = eval.bspline%*%beta_coef
 
beta = matrix(ncol=length(beta_line)*2, nrow=length(beta_line))
for(i in 1:length(beta_line)) beta[i,] = c(rep(0, (i-1)), beta_line, rep(0,length(beta_line)-i+1))
beta = beta[,(length(beta_line)+1):(length(beta_line)*2)]

if(plots){
  matplot(cbind(delta1,alpha), type="l")
  matplot(cbind(delta2,gamma,beta_line), type="l")
  image(t(beta), col  = gray((0:32)/32))
  }

alphabeta = matrix(alpha, byrow=T, ncol=length(alpha), nrow=length(alpha))*beta # In beta, the columns are 's', while the rows are 't'
alphabeta_integral  = rowSums(alphabeta)*(timevec[2]-timevec[1])  # Integral of ab-function: ab(t) = \int af(s)bf(t,s)ds dt gives ab-path
```

Setting seed:

```r
if(cluster) RNGkind("L'Ecuyer-CMRG")
seed = 340953
set.seed(seed)
```


### Simulation 4
In this simulation, $M$ depends on $Z$ and $X$, $Y$ only depends on $M$ and not on $X$ (there is complete mediation).


```r
Pall  = matrix(0, iters, (2*len^2+5*len)) #2 surfaces (abs, bf), 5 functions (abf,af,gf,d1,d2)
all   = matrix(0, iters, (2*len^2+5*len))

for(i in 1:iters){
  X          = rbinom(sub,1,0.5)
  
  epsilon    = matrix(rnorm(sub*len,0,sd_epsilon), ncol=len,nrow=(sub))

  design_M = cbind(int=1,X=X)
  coeff_M  = cbind(delta1,alpha)
  colnames(coeff_M) = c("delta1","alpha")

  M    = t(tcrossprod(coeff_M, design_M))+epsilon

  eta  = matrix(rnorm(sub*len,0,sd_eta), ncol=len,nrow=(sub))

  Y    = t(sapply(1:sub, function(i) create_functional_y(x=X[i], m=M[i,], etai=eta[i,], delta2=delta2, gamma=gamma, beta=beta)))

  dta  = data.frame(X=X,M=M,Y=Y)
  
  out <- sff_Mediation(X,Y,t(M),mediatorMethod="fosr2s", nbasis=nbasis,norder=norder,lambda=lambda, plot=plots, boot=TRUE)
  
  if(cluster==TRUE){
    time1 = Sys.time()
    stat = boot(data = dta, statistic = fbootstrap_ML,R = nboot,sim = "ordinary",stype = "i",parallel = "multicore",ncpus=15)
    time2 = Sys.time()
    }else{
      time1 = Sys.time()
      stat = boot(data = dta, statistic = fbootstrap_ML,R = nboot,sim = "ordinary",stype = "i")
      time2 = Sys.time()
    }

  print(time2-time1)
  
  # Calculating p-value
  Pall[i,] = 2*apply(cbind(colMeans(stat$t < 2*matrix(stat$t0, ncol=ncol(stat$t), nrow=nrow(stat$t), byrow = TRUE), na.rm = TRUE),colMeans(stat$t > 2*matrix(stat$t0, ncol=ncol(stat$t), nrow=nrow(stat$t), byrow = TRUE), na.rm = TRUE)),1, min)
  all[i,] = stat$t0
  }
```

```
## Time difference of 2.118572 mins
## Time difference of 2.846505 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5, 14, 4, 9, 12, 10 encountered errors in user code, all
## values of the jobs will be affected
```

```
## Time difference of 2.05405 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 12, 6, 15, 5, 3, 13, 9, 10 encountered errors in user
## code, all values of the jobs will be affected
```

```
## Time difference of 1.882716 mins
## Time difference of 2.660101 mins
## Time difference of 3.122904 mins
## Time difference of 3.607651 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1, 10 encountered error in user code, all values of the job
## will be affected
```

```
## Time difference of 1.869338 mins
## Time difference of 2.496209 mins
## Time difference of 1.989512 mins
## Time difference of 3.899521 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 12 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Time difference of 2.02388 mins
## Time difference of 2.843786 mins
## Time difference of 2.529254 mins
## Time difference of 3.167641 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1, 5 encountered error in user code, all values of the job
## will be affected
```

```
## Time difference of 3.603618 mins
## Time difference of 1.997791 mins
## Time difference of 2.652989 mins
## Time difference of 1.8417 mins
## Time difference of 1.916056 mins
## Time difference of 1.921589 mins
## Time difference of 2.259942 mins
## Time difference of 1.818497 mins
## Time difference of 1.828306 mins
## Time difference of 1.887138 mins
## Time difference of 1.937127 mins
## Time difference of 2.288485 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 8, 7 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Time difference of 2.406618 mins
## Time difference of 3.395087 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 15 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Time difference of 2.319243 mins
## Time difference of 2.242903 mins
## Time difference of 2.897562 mins
## Time difference of 3.510803 mins
## Time difference of 2.165662 mins
## Time difference of 3.73203 mins
## Time difference of 3.549927 mins
## Time difference of 3.090297 mins
## Time difference of 2.966601 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 10, 14 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Time difference of 3.763608 mins
## Time difference of 2.887747 mins
## Time difference of 2.29289 mins
## Time difference of 1.905196 mins
## Time difference of 2.408077 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 15 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Time difference of 3.235079 mins
## Time difference of 3.193563 mins
## Time difference of 2.406425 mins
## Time difference of 2.25626 mins
## Time difference of 2.359931 mins
## Time difference of 2.788895 mins
## Time difference of 2.523934 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5, 4, 6, 9, 8, 15, 7, 11 encountered errors in user code,
## all values of the jobs will be affected
```

```
## Time difference of 1.901564 mins
## Time difference of 2.679906 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 7 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Time difference of 3.650852 mins
## Time difference of 4.686613 mins
## Time difference of 2.598039 mins
## Time difference of 2.367073 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 10, 6, 15, 2, 3, 1, 7, 8 encountered errors in user code,
## all values of the jobs will be affected
```

```
## Time difference of 1.763932 mins
## Time difference of 1.809773 mins
## Time difference of 2.907531 mins
## Time difference of 2.259825 mins
## Time difference of 2.151958 mins
## Time difference of 1.868728 mins
## Time difference of 3.321926 mins
## Time difference of 2.566691 mins
## Time difference of 2.581793 mins
## Time difference of 2.121748 mins
## Time difference of 2.693639 mins
## Time difference of 1.821573 mins
## Time difference of 2.298634 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 13 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Time difference of 2.576946 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Time difference of 2.766055 mins
## Time difference of 2.503187 mins
## Time difference of 2.334478 mins
## Time difference of 2.254184 mins
## Time difference of 2.083178 mins
## Time difference of 2.109641 mins
## Time difference of 3.542474 mins
## Time difference of 2.33703 mins
## Time difference of 2.158453 mins
## Time difference of 2.795425 mins
## Time difference of 3.549499 mins
## Time difference of 1.825318 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 7, 4, 15 encountered errors in user code, all values of
## the jobs will be affected
```

```
## Time difference of 2.267138 mins
## Time difference of 3.092901 mins
## Time difference of 2.076754 mins
## Time difference of 2.540567 mins
## Time difference of 3.347799 mins
## Time difference of 2.364227 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 13 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Time difference of 3.500478 mins
## Time difference of 2.443364 mins
## Time difference of 2.356214 mins
## Time difference of 2.466675 mins
## Time difference of 2.348205 mins
## Time difference of 1.987015 mins
## Time difference of 3.251234 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 3 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Time difference of 2.964665 mins
## Time difference of 2.684545 mins
## Time difference of 2.776675 mins
## Time difference of 2.888261 mins
## Time difference of 2.22037 mins
## Time difference of 2.404761 mins
## Time difference of 1.937743 mins
## Time difference of 2.662264 mins
## Time difference of 2.371066 mins
## Time difference of 1.811546 mins
## Time difference of 3.522473 mins
## Time difference of 2.519863 mins
## Time difference of 4.479725 mins
## Time difference of 2.046017 mins
## Time difference of 2.890974 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 3 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Time difference of 3.913011 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 3 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Time difference of 3.157695 mins
## Time difference of 2.146695 mins
## Time difference of 3.726693 mins
## Time difference of 1.828857 mins
## Time difference of 1.530557 mins
## Time difference of 3.226795 mins
## Time difference of 2.990816 mins
## Time difference of 1.920709 mins
## Time difference of 1.919698 mins
## Time difference of 2.089314 mins
## Time difference of 2.616557 mins
## Time difference of 2.41229 mins
## Time difference of 2.037893 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5, 11, 14, 8, 13, 9, 12, 6, 10, 4, 3, 1 encountered errors
## in user code, all values of the jobs will be affected
```

```
## Time difference of 1.239029 mins
## Time difference of 3.514501 mins
## Time difference of 3.232088 mins
## Time difference of 3.245278 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 8 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Time difference of 3.955404 mins
## Time difference of 3.719552 mins
## Time difference of 3.022316 mins
## Time difference of 2.355876 mins
## Time difference of 1.975129 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 7 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Time difference of 1.998927 mins
## Time difference of 2.393319 mins
## Time difference of 3.979538 mins
## Time difference of 2.38099 mins
## Time difference of 2.299805 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Time difference of 2.208552 mins
## Time difference of 2.201277 mins
## Time difference of 3.305351 mins
## Time difference of 2.321615 mins
## Time difference of 1.933833 mins
## Time difference of 2.087338 mins
## Time difference of 2.321096 mins
## Time difference of 2.727454 mins
## Time difference of 2.376436 mins
## Time difference of 1.967285 mins
## Time difference of 3.715027 mins
## Time difference of 3.391173 mins
## Time difference of 2.499426 mins
## Time difference of 3.678502 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 10, 7 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Time difference of 1.798536 mins
## Time difference of 2.97865 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus): all
## scheduled cores encountered errors in user code
```

```
## Time difference of 1.104217 mins
## Time difference of 1.885146 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 14 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Time difference of 2.726762 mins
## Time difference of 3.70984 mins
## Time difference of 2.322952 mins
## Time difference of 2.630645 mins
## Time difference of 2.147395 mins
## Time difference of 3.921089 mins
## Time difference of 3.019991 mins
## Time difference of 2.506857 mins
## Time difference of 2.555701 mins
## Time difference of 1.71345 mins
## Time difference of 1.89431 mins
## Time difference of 1.977102 mins
## Time difference of 2.630853 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 8, 9, 11, 4 encountered errors in user code, all values of
## the jobs will be affected
```

```
## Time difference of 3.474157 mins
## Time difference of 2.860514 mins
## Time difference of 2.711687 mins
## Time difference of 3.43528 mins
## Time difference of 2.72639 mins
## Time difference of 2.507837 mins
## Time difference of 1.885436 mins
## Time difference of 2.932612 mins
## Time difference of 3.244894 mins
## Time difference of 2.656083 mins
## Time difference of 2.419581 mins
## Time difference of 1.98335 mins
## Time difference of 3.62975 mins
## Time difference of 3.010969 mins
## Time difference of 2.038008 mins
## Time difference of 3.455686 mins
## Time difference of 2.153796 mins
## Time difference of 2.479413 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 14, 5, 3, 1, 11 encountered errors in user code, all
## values of the jobs will be affected
```

```
## Time difference of 3.437506 mins
## Time difference of 2.79527 mins
## Time difference of 3.184222 mins
## Time difference of 3.053546 mins
## Time difference of 2.524819 mins
## Time difference of 2.306792 mins
## Time difference of 2.399605 mins
## Time difference of 3.69194 mins
## Time difference of 2.243736 mins
## Time difference of 2.286006 mins
## Time difference of 3.924815 mins
## Time difference of 2.438938 mins
## Time difference of 2.119596 mins
## Time difference of 2.859286 mins
## Time difference of 2.441379 mins
## Time difference of 1.792122 mins
## Time difference of 1.664231 mins
## Time difference of 2.136459 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 10, 14, 15, 7, 9, 12 encountered errors in user code, all
## values of the jobs will be affected
```

```
## Time difference of 1.511959 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 13 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Time difference of 1.736775 mins
## Time difference of 2.983372 mins
## Time difference of 1.693086 mins
## Time difference of 3.047126 mins
## Time difference of 2.312614 mins
## Time difference of 1.82226 mins
## Time difference of 2.088169 mins
## Time difference of 2.673966 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 12, 5, 6 encountered errors in user code, all values of
## the jobs will be affected
```

```
## Time difference of 4.103208 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5, 4, 10, 14, 1, 13, 8, 11 encountered errors in user
## code, all values of the jobs will be affected
```

```
## Time difference of 2.013921 mins
## Time difference of 2.251634 mins
## Time difference of 3.568339 mins
## Time difference of 3.851445 mins
## Time difference of 2.771608 mins
## Time difference of 1.966095 mins
## Time difference of 2.659072 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 6 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Time difference of 1.686495 mins
## Time difference of 2.250239 mins
## Time difference of 1.722229 mins
## Time difference of 3.394157 mins
## Time difference of 3.897701 mins
## Time difference of 1.798585 mins
## Time difference of 2.16062 mins
## Time difference of 2.286757 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 8 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Time difference of 2.215561 mins
## Time difference of 2.121416 mins
## Time difference of 2.591542 mins
## Time difference of 1.8276 mins
## Time difference of 3.045507 mins
## Time difference of 2.126461 mins
## Time difference of 1.660271 mins
## Time difference of 2.018989 mins
## Time difference of 2.670879 mins
## Time difference of 2.536754 mins
## Time difference of 2.524921 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 15 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Time difference of 2.164734 mins
## Time difference of 2.524444 mins
## Time difference of 3.009522 mins
## Time difference of 2.986011 mins
## Time difference of 2.285432 mins
## Time difference of 1.911118 mins
## Time difference of 2.455779 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5, 15, 1 encountered errors in user code, all values of
## the jobs will be affected
```

```
## Time difference of 2.767709 mins
## Time difference of 4.723745 mins
## Time difference of 2.170259 mins
## Time difference of 3.105345 mins
## Time difference of 2.781131 mins
## Time difference of 3.350401 mins
## Time difference of 1.771752 mins
## Time difference of 2.121992 mins
## Time difference of 2.668265 mins
## Time difference of 2.51968 mins
## Time difference of 3.12202 mins
## Time difference of 3.406642 mins
## Time difference of 1.873358 mins
## Time difference of 2.620569 mins
## Time difference of 1.727988 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9, 1, 4, 2, 8, 3 encountered errors in user code, all
## values of the jobs will be affected
```

```
## Time difference of 2.796353 mins
## Time difference of 2.091699 mins
## Time difference of 2.387883 mins
## Time difference of 2.445339 mins
## Time difference of 1.741361 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 7 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Time difference of 3.881957 mins
## Time difference of 3.025668 mins
## Time difference of 2.886898 mins
## Time difference of 2.444635 mins
## Time difference of 2.087982 mins
## Time difference of 1.870962 mins
## Time difference of 2.10776 mins
## Time difference of 3.364175 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 14 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Time difference of 3.644575 mins
## Time difference of 2.848903 mins
## Time difference of 2.366343 mins
## Time difference of 2.398451 mins
## Time difference of 3.261322 mins
## Time difference of 2.797816 mins
## Time difference of 2.987978 mins
## Time difference of 2.842085 mins
## Time difference of 2.60595 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 14 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Time difference of 2.563942 mins
## Time difference of 2.975711 mins
## Time difference of 3.087595 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9, 4, 12, 11, 13 encountered errors in user code, all
## values of the jobs will be affected
```

```
## Time difference of 2.000108 mins
## Time difference of 3.531183 mins
## Time difference of 2.884275 mins
## Time difference of 2.552172 mins
## Time difference of 2.911324 mins
## Time difference of 2.24467 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 14, 4, 11, 8, 3, 13, 9, 15, 12, 5, 7 encountered errors in
## user code, all values of the jobs will be affected
```

```
## Time difference of 3.083311 mins
## Time difference of 1.99837 mins
## Time difference of 2.390493 mins
## Time difference of 2.659622 mins
## Time difference of 4.1123 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 12 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Time difference of 2.796675 mins
## Time difference of 4.064622 mins
## Time difference of 2.69671 mins
## Time difference of 2.778639 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2, 14, 10, 13, 5, 8, 12, 1, 3, 15, 4, 11 encountered
## errors in user code, all values of the jobs will be affected
```

```
## Time difference of 2.207825 mins
## Time difference of 3.205934 mins
## Time difference of 2.256446 mins
## Time difference of 2.506141 mins
## Time difference of 2.771161 mins
## Time difference of 2.303625 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 8, 2, 7, 5, 6, 12, 1, 3 encountered errors in user code,
## all values of the jobs will be affected
```

```
## Time difference of 1.558145 mins
## Time difference of 2.311811 mins
## Time difference of 2.266165 mins
## Time difference of 2.018498 mins
## Time difference of 2.840997 mins
## Time difference of 2.999306 mins
## Time difference of 3.778727 mins
## Time difference of 2.435546 mins
## Time difference of 2.326718 mins
## Time difference of 3.539915 mins
## Time difference of 2.236052 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Time difference of 5.097709 mins
## Time difference of 2.335472 mins
## Time difference of 2.321727 mins
## Time difference of 1.723848 mins
## Time difference of 3.820813 mins
## Time difference of 3.246968 mins
## Time difference of 2.174701 mins
## Time difference of 2.242181 mins
## Time difference of 1.960498 mins
## Time difference of 1.726679 mins
## Time difference of 2.408706 mins
## Time difference of 3.381123 mins
## Time difference of 3.461523 mins
## Time difference of 2.837565 mins
## Time difference of 3.022659 mins
## Time difference of 3.319508 mins
## Time difference of 3.957159 mins
## Time difference of 2.517341 mins
## Time difference of 4.420325 mins
## Time difference of 2.618607 mins
## Time difference of 2.505783 mins
## Time difference of 3.270861 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Time difference of 2.281684 mins
## Time difference of 3.140785 mins
## Time difference of 2.743359 mins
## Time difference of 2.272373 mins
## Time difference of 4.45897 mins
## Time difference of 2.588781 mins
## Time difference of 2.247027 mins
## Time difference of 3.558531 mins
## Time difference of 2.034659 mins
## Time difference of 1.818677 mins
## Time difference of 2.384067 mins
## Time difference of 2.712624 mins
## Time difference of 3.297997 mins
## Time difference of 3.958129 mins
## Time difference of 2.059059 mins
## Time difference of 3.2276 mins
## Time difference of 1.946797 mins
## Time difference of 3.340076 mins
## Time difference of 1.835995 mins
## Time difference of 2.982347 mins
## Time difference of 2.098402 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 3 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Time difference of 2.287799 mins
## Time difference of 2.468379 mins
## Time difference of 2.15633 mins
## Time difference of 3.106903 mins
## Time difference of 2.514578 mins
## Time difference of 4.507681 mins
## Time difference of 2.410772 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2, 10, 13 encountered errors in user code, all values of
## the jobs will be affected
```

```
## Time difference of 2.691764 mins
## Time difference of 3.204015 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Time difference of 1.854721 mins
## Time difference of 3.461064 mins
## Time difference of 3.317348 mins
## Time difference of 2.609455 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Time difference of 2.513584 mins
## Time difference of 2.285474 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 6 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Time difference of 2.946353 mins
## Time difference of 3.736854 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 8 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Time difference of 2.113852 mins
## Time difference of 2.855117 mins
## Time difference of 3.592899 mins
## Time difference of 2.460526 mins
## Time difference of 2.029267 mins
## Time difference of 2.783168 mins
## Time difference of 1.982701 mins
## Time difference of 3.030247 mins
## Time difference of 5.228361 mins
## Time difference of 2.902945 mins
## Time difference of 2.59712 mins
## Time difference of 2.972083 mins
## Time difference of 2.08532 mins
## Time difference of 2.641071 mins
## Time difference of 2.062706 mins
## Time difference of 2.456977 mins
## Time difference of 4.191547 mins
## Time difference of 3.731623 mins
## Time difference of 2.032369 mins
## Time difference of 2.743915 mins
## Time difference of 2.131892 mins
## Time difference of 2.705711 mins
## Time difference of 2.336005 mins
## Time difference of 1.820571 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5, 15, 11, 3, 12, 14, 9, 13, 4, 7, 10, 2, 8, 1 encountered
## errors in user code, all values of the jobs will be affected
```

```
## Time difference of 1.66297 mins
## Time difference of 2.591705 mins
## Time difference of 2.489133 mins
## Time difference of 2.932688 mins
## Time difference of 2.906796 mins
## Time difference of 1.96501 mins
## Time difference of 4.524434 mins
## Time difference of 1.907593 mins
## Time difference of 3.211259 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Time difference of 1.703178 mins
## Time difference of 3.455115 mins
## Time difference of 3.016233 mins
## Time difference of 2.228012 mins
## Time difference of 2.289682 mins
## Time difference of 2.403241 mins
## Time difference of 2.416036 mins
## Time difference of 2.937152 mins
## Time difference of 2.068133 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 10, 2, 1, 12, 5, 14, 8, 9, 7, 6 encountered errors in user
## code, all values of the jobs will be affected
```

```
## Time difference of 2.483582 mins
## Time difference of 1.801003 mins
## Time difference of 2.301076 mins
## Time difference of 2.783145 mins
## Time difference of 2.100197 mins
## Time difference of 1.869497 mins
## Time difference of 2.833315 mins
## Time difference of 2.938002 mins
## Time difference of 3.04849 mins
## Time difference of 1.943277 mins
## Time difference of 1.89194 mins
## Time difference of 3.086117 mins
## Time difference of 1.85883 mins
## Time difference of 2.6166 mins
## Time difference of 2.167113 mins
## Time difference of 2.488259 mins
## Time difference of 3.868209 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 8 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Time difference of 2.397979 mins
## Time difference of 2.051089 mins
## Time difference of 1.797232 mins
## Time difference of 2.153121 mins
## Time difference of 2.150155 mins
## Time difference of 1.822457 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 12 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Time difference of 1.937265 mins
## Time difference of 3.23046 mins
## Time difference of 2.842138 mins
## Time difference of 2.59229 mins
## Time difference of 2.532604 mins
## Time difference of 1.789318 mins
## Time difference of 2.372344 mins
## Time difference of 3.208592 mins
## Time difference of 2.632266 mins
## Time difference of 3.211841 mins
## Time difference of 3.432707 mins
## Time difference of 2.266558 mins
## Time difference of 3.700309 mins
## Time difference of 2.501621 mins
## Time difference of 2.919656 mins
## Time difference of 2.957065 mins
## Time difference of 2.000744 mins
## Time difference of 2.537742 mins
## Time difference of 3.108143 mins
## Time difference of 2.688378 mins
## Time difference of 2.527475 mins
## Time difference of 3.901713 mins
## Time difference of 2.693378 mins
## Time difference of 2.01823 mins
## Time difference of 2.078798 mins
## Time difference of 2.012852 mins
## Time difference of 2.799076 mins
## Time difference of 2.380365 mins
## Time difference of 1.864487 mins
## Time difference of 2.185915 mins
## Time difference of 1.913002 mins
## Time difference of 1.999331 mins
## Time difference of 1.540279 mins
## Time difference of 2.944531 mins
## Time difference of 2.808069 mins
## Time difference of 2.483009 mins
## Time difference of 2.211458 mins
## Time difference of 2.771363 mins
## Time difference of 2.283705 mins
## Time difference of 1.89826 mins
## Time difference of 2.649281 mins
## Time difference of 3.02878 mins
## Time difference of 2.335477 mins
## Time difference of 2.342212 mins
## Time difference of 3.629857 mins
## Time difference of 3.604345 mins
## Time difference of 3.408809 mins
## Time difference of 3.670494 mins
## Time difference of 2.338006 mins
## Time difference of 3.263768 mins
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 13, 7, 4, 2, 5, 3, 12 encountered errors in user code, all
## values of the jobs will be affected
```

```
## Time difference of 4.053519 mins
## Time difference of 3.076918 mins
## Time difference of 2.561433 mins
## Time difference of 2.682963 mins
## Time difference of 1.846649 mins
```

```r
colnames(Pall) = colnames(all) = colnames(stat$t) = names(stat$t0)

if(plots){
  sff_vec2list(Pall[1,], stats="bsurface", plot=TRUE, returns=FALSE)
  matplot(sapply(1:10, function(i) sff_vec2list(stat$t[i,], stats="abfunction")), type="l")
}
save(all,Pall, file=paste0("SFF_Simulation_4_",format(Sys.time(), "%Y%m%d-%H%M"),".Rdata"))
```


### Null Simulation 1. 
In this simulation, $M$ is pure noise (and is not dependent on $X$), and $Y$ only depends on $X$.


```r
Pall  = matrix(0, iters, (2*len^2+5*len)) #2 surfaces (abs, bf), 5 functions (abf,af,gf,d1,d2)
all   = matrix(0, iters, (2*len^2+5*len))

for(i in 1:iters){
  X          = rbinom(sub,1,0.5)
  
  epsilon    = matrix(rnorm(sub*len,0,sd_epsilon), ncol=len,nrow=(sub))

  M    = epsilon

  eta  = matrix(rnorm(sub*len,0,sd_eta), ncol=len,nrow=(sub))

  Y    = t(sapply(1:sub, function(i) create_functional_y(x=X[i], m=M[i,], etai=eta[i,], delta2=delta2, gamma=gamma, beta=matrix(0, ncol=ncol(beta), nrow=nrow(beta)))))

 dta  = data.frame(X=X,M=M,Y=Y)
  
  out  = sff_Mediation(X,Y,t(M),mediatorMethod="fosr2s", nbasis=nbasis,norder=norder,lambda=lambda, plot=plots, boot=TRUE)
  
  if(cluster==TRUE){
    time1 = Sys.time()
    stat = boot(data = dta, statistic = fbootstrap_ML,R = nboot,sim = "ordinary",stype = "i",parallel = "multicore",ncpus=15)
    time2 = Sys.time()
    }else{
    stat = boot(data = dta, statistic = fbootstrap_ML,R = nboot,sim = "ordinary",stype = "i")
    }
                
  
  # Calculating p-value
  Pall[i,] = 2*apply(cbind(colMeans(stat$t < 2*matrix(stat$t0, ncol=ncol(stat$t), nrow=nrow(stat$t), byrow = TRUE), na.rm = TRUE),colMeans(stat$t > 2*matrix(stat$t0, ncol=ncol(stat$t), nrow=nrow(stat$t), byrow = TRUE), na.rm = TRUE)),1, min)
  all[i,] = stat$t0
  }
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 11, 7 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 15, 13, 1, 5 encountered errors in user code, all values
## of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 10, 11, 5 encountered errors in user code, all values of
## the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1, 2, 4, 15, 11, 6, 14, 9, 5 encountered error in user
## code, all values of the job will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9, 13, 4 encountered errors in user code, all values of
## the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 3, 4, 9 encountered errors in user code, all values of the
## jobs will be affected
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
## scheduled core 1 encountered error in user code, all values of the job
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 7, 1, 15 encountered errors in user code, all values of
## the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4, 5, 9, 6 encountered errors in user code, all values of
## the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 13, 10 encountered errors in user code, all values of the
## jobs will be affected
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
## scheduled cores 14, 10 encountered errors in user code, all values of the
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
## scheduled cores 3, 8, 2 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 14, 12, 1 encountered errors in user code, all values of
## the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 8, 9, 14, 1, 5, 7, 6 encountered errors in user code, all
## values of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 3, 13, 12, 5, 2 encountered errors in user code, all
## values of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 10, 14 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 15 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4, 1 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 7, 12 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 14 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5, 2, 14, 12 encountered errors in user code, all values
## of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 13, 12, 14 encountered errors in user code, all values of
## the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2, 1, 11, 13, 3, 7, 12, 9, 4 encountered errors in user
## code, all values of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 10 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9, 11, 12, 14, 1 encountered errors in user code, all
## values of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1, 9, 4, 13, 12 encountered error in user code, all values
## of the job will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5, 6, 4 encountered errors in user code, all values of the
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
## scheduled cores 3 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 14, 10 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1, 14 encountered error in user code, all values of the job
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 15, 1 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 13, 1, 5, 3, 14, 11 encountered errors in user code, all
## values of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 8, 11, 3 encountered errors in user code, all values of
## the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus): all
## scheduled cores encountered errors in user code
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 3, 7, 1, 11, 8, 9, 10, 13, 12, 15 encountered errors in
## user code, all values of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 12 encountered errors in user code, all values of the jobs
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

```r
colnames(Pall) = colnames(all) = colnames(stat$t) = names(stat$t0)

if(plots){
  sff_vec2list(Pall[1,], stats="bsurface", plot=TRUE, returns=FALSE)
  matplot(sapply(1:10, function(i) sff_vec2list(stat$t[i,], stats="abfunction")), type="l")
}

save(all,Pall, file=paste0("SFF_Simulation_1_",format(Sys.time(), "%Y%m%d-%H%M"),".Rdata"))
```



### Null Simulation 2. 
In this simulation, $M$ depends on $X$, but $Y$ only depends on $X$ and not on $M$.


```r
Pall  = matrix(0, iters, (2*len^2+5*len)) #2 surfaces (abs, bf), 5 functions (abf,af,gf,d1,d2)
all   = matrix(0, iters, (2*len^2+5*len))

for(i in 1:iters){
  X          = rbinom(sub,1,0.5)
  
  epsilon    = matrix(rnorm(sub*len,0,sd_epsilon), ncol=len,nrow=(sub))

  design_M = cbind(int=1,X=X)
  coeff_M  = cbind(delta1,alpha)
  colnames(coeff_M) = c("delta1","alpha")

  M    = t(tcrossprod(coeff_M, design_M))+epsilon

  eta  = matrix(rnorm(sub*len,0,sd_eta), ncol=len,nrow=(sub))

  Y    = t(sapply(1:sub, function(i) create_functional_y(x=X[i], m=M[i,], etai=eta[i,], delta2=delta2, gamma=gamma, beta=matrix(0, ncol=ncol(beta), nrow=nrow(beta)))))

  dta  = data.frame(X=X,M=M,Y=Y)
 
  out <- sff_Mediation(X,Y,t(M),mediatorMethod="fosr2s", nbasis=nbasis,norder=norder,lambda=lambda, plot=plots, boot=TRUE) 

  if(cluster==TRUE){
    time1 = Sys.time()
    stat = boot(data = dta, statistic = fbootstrap_ML,R = nboot,sim = "ordinary",stype = "i",parallel = "multicore",ncpus=15)
    time2 = Sys.time()
    }else{
    stat = boot(data = dta, statistic = fbootstrap_ML,R = nboot,sim = "ordinary",stype = "i")
    }

                
  
  # Calculating p-value
  Pall[i,] = 2*apply(cbind(colMeans(stat$t < 2*matrix(stat$t0, ncol=ncol(stat$t), nrow=nrow(stat$t), byrow = TRUE), na.rm = TRUE),colMeans(stat$t > 2*matrix(stat$t0, ncol=ncol(stat$t), nrow=nrow(stat$t), byrow = TRUE), na.rm = TRUE)),1, min)
  all[i,] = stat$t0
  }
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 3, 5, 11, 7 encountered errors in user code, all values of
## the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 11, 14, 2 encountered errors in user code, all values of
## the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 3 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 12 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 13, 14, 9 encountered errors in user code, all values of
## the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9, 11, 2, 1 encountered errors in user code, all values of
## the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 13, 7, 5, 6, 12 encountered errors in user code, all
## values of the jobs will be affected
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
## scheduled cores 2, 9 encountered errors in user code, all values of the
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
## scheduled cores 11 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 13, 4 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 7 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 14, 8 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 6, 7, 10, 11, 3, 1, 12 encountered errors in user code,
## all values of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5, 1 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1, 5 encountered error in user code, all values of the job
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 3 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1, 4, 2, 11, 15 encountered error in user code, all values
## of the job will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5, 1, 12 encountered errors in user code, all values of
## the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 12 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4, 9 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 6 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 7, 15, 6 encountered errors in user code, all values of
## the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1, 15 encountered error in user code, all values of the job
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 12, 3, 6 encountered errors in user code, all values of
## the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 6 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 14 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 14, 13, 3 encountered errors in user code, all values of
## the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 11, 13, 10, 3, 8, 1, 5, 15 encountered errors in user
## code, all values of the jobs will be affected
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
## scheduled cores 15, 2, 7, 4 encountered errors in user code, all values of
## the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1 encountered error in user code, all values of the job
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2, 14, 12 encountered errors in user code, all values of
## the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9, 5, 12, 4, 1, 6, 13, 7 encountered errors in user code,
## all values of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 3, 11 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 11, 3 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1, 6 encountered error in user code, all values of the job
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 11 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1 encountered error in user code, all values of the job
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 15, 12, 5 encountered errors in user code, all values of
## the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 15 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 14 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 3, 8, 11, 5, 10 encountered errors in user code, all
## values of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 10, 14, 15 encountered errors in user code, all values of
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
## scheduled cores 14 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 13, 14, 5 encountered errors in user code, all values of
## the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1, 10 encountered error in user code, all values of the job
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5, 3 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 12 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 7, 13, 11 encountered errors in user code, all values of
## the jobs will be affected
```

```r
colnames(Pall) = colnames(all) = colnames(stat$t) = names(stat$t0)

if(plots){
  sff_vec2list(Pall[1,], stats="bsurface", plot=TRUE, returns=FALSE)
  matplot(sapply(1:10, function(i) sff_vec2list(stat$t[i,], stats="abfunction")), type="l")
}

save(all,Pall, file=paste0("SFF_Simulation_2_",format(Sys.time(), "%Y%m%d-%H%M"),".Rdata"))
```


### Null simulation 3
In this simulation, $M$ depends on $\delta_1$ for all subjects (and does not depend on $X$), $Y$ only depends on $M$ and not on $X$.


```r
Pall  = matrix(0, iters, (2*len^2+5*len)) #2 surfaces (abs, bf), 5 functions (abf,af,gf,d1,d2)
all   = matrix(0, iters, (2*len^2+5*len))

for(i in 1:iters){
  X          = rbinom(sub,1,0.5)
  
  epsilon    = matrix(rnorm(sub*len,0,sd_epsilon), ncol=len,nrow=(sub))

  design_M = rep(1,sub)
  coeff_M  = cbind(delta1)
  colnames(coeff_M) = c("delta1")

  M    = t(tcrossprod(coeff_M, design_M))+epsilon

  eta  = matrix(rnorm(sub*len,0,sd_eta), ncol=len,nrow=(sub))

  Y    = t(sapply(1:sub, function(i) create_functional_y(x=X[i], m=M[i,], etai=eta[i,], delta2=delta2, gamma=rep(0, length(gamma)), beta=beta)))

  dta  = data.frame(X=X,M=M,Y=Y)
  
  out <- sff_Mediation(X,Y,t(M),mediatorMethod="fosr2s",nbasis=nbasis,norder=norder,lambda=lambda, plot=plots, boot=TRUE)

  if(cluster==TRUE){
    time1 = Sys.time()
    stat = boot(data = dta, statistic = fbootstrap_ML,R = nboot,sim = "ordinary",stype = "i",parallel = "multicore",ncpus=15)
    time2 = Sys.time()
    }else{
    stat = boot(data = dta, statistic = fbootstrap_ML,R = nboot,sim = "ordinary",stype = "i")
    }

                
  
  # Calculating p-value
  Pall[i,] = 2*apply(cbind(colMeans(stat$t < 2*matrix(stat$t0, ncol=ncol(stat$t), nrow=nrow(stat$t), byrow = TRUE), na.rm = TRUE),colMeans(stat$t > 2*matrix(stat$t0, ncol=ncol(stat$t), nrow=nrow(stat$t), byrow = TRUE), na.rm = TRUE)),1, min)
  all[i,] = stat$t0
  }
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 11 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9, 1, 14, 6, 4, 2, 12, 11, 5 encountered errors in user
## code, all values of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9, 2, 15, 10, 11, 3, 13, 14, 8 encountered errors in user
## code, all values of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4, 6 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 4 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 12 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 13, 6, 9, 1 encountered errors in user code, all values of
## the jobs will be affected
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
## scheduled cores 2, 3 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 2, 15 encountered errors in user code, all values of the
## jobs will be affected
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
## scheduled cores 4, 6 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus): all
## scheduled cores encountered errors in user code
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 3 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 12, 10, 2, 14, 7 encountered errors in user code, all
## values of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 7 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 11, 14 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1 encountered error in user code, all values of the job
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 15 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled core 1, 15, 11, 5 encountered error in user code, all values of
## the job will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5, 3, 8, 12 encountered errors in user code, all values of
## the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 13 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 15 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 10, 11 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 12, 11, 14 encountered errors in user code, all values of
## the jobs will be affected
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
## scheduled core 1, 5, 4, 11, 14, 7, 13, 6 encountered error in user code,
## all values of the job will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 5 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 7 encountered errors in user code, all values of the jobs
## will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9, 6, 7, 5, 11, 3, 10, 15, 2, 8, 13 encountered errors in
## user code, all values of the jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 9, 1 encountered errors in user code, all values of the
## jobs will be affected
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
## scheduled cores 6, 7, 4 encountered errors in user code, all values of the
## jobs will be affected
```

```
## Warning in parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus):
## scheduled cores 7, 2, 3, 11, 1, 13, 5 encountered errors in user code, all
## values of the jobs will be affected
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

```r
colnames(Pall) = colnames(all) = colnames(stat$t) = names(stat$t0)

if(plots){
  sff_vec2list(Pall[1,], stats="bsurface", plot=TRUE, returns=FALSE)
  matplot(sapply(1:10, function(i) sff_vec2list(stat$t[i,], stats="abfunction")), type="l")
}

save(all,Pall, file=paste0("SFF_Simulation_3_",format(Sys.time(), "%Y%m%d-%H%M"),".Rdata"))
```
