---
title: "Data Application in Cluster"
author: "Yenny Webb-Vargas"
date: "Wednesday, February 04, 2015"
output: html_document
---


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
  meta_file = '/home/bst/student/ywebbvar/Github/functional_mediation/data_application/meta.mat'
  }else{
  source('~/GitHub/functional_mediation/sfs_Mediation.R')  
  source('~/GitHub/functional_mediation/est_se_fgam.R')  
  meta_file = '~/GitHub/functional_mediation/data_application/meta.mat'
  }
```



```r
library(R.matlab)
```

```
## R.matlab v3.1.1 (2014-10-10) successfully loaded. See ?R.matlab for help.
## 
## Attaching package: 'R.matlab'
## 
## The following objects are masked from 'package:base':
## 
##     getOption, isOpen
```

```r
meta = readMat(meta_file, fixNames=TRUE)[['meta']]
```


```r
list_names = 'heat_crating: {1x93 cell}
         soundpain_crating: {1x93 cell}
              iads_crating: {1x93 cell}
          heat_orating_int: {1x93 cell}
          heat_orating_unp: {1x93 cell}
     soundpain_orating_int: {1x93 cell}
     sounpdain_orating_unp: {1x93 cell}
          iads_orating_int: {1x93 cell}
          iads_orating_unp: {1x93 cell}
                       run: {1x93 cell}
                     trial: {1x93 cell}
                 condition: {1x93 cell}
                  iads_set: {1x93 cell}
                heat_level: {1x93 cell}
         heat_crating_peak: {1x93 cell}
           soundpain_level: {1x93 cell}
    soundpain_crating_peak: {1x93 cell}
                iads_level: {1x93 cell}
         iads_crating_peak: {1x93 cell}
                  subjects: {1x93 cell}
                   ratings: {1x93 cell}
                      temp: {1x93 cell}'

list_names = gsub(' {1x93 cell}\n', '', list_names, fixed=TRUE)
list_names = gsub(' {1x93 cell}', '', list_names, fixed=TRUE)
list_names = strsplit(list_names, ":")[[1]]
list_names = sapply(list_names, function(x) gsub(' ', '', x))

names(meta) = list_names
```

Dropping all cases with missing mediator and/or missing outcome:

```r
get_data <- function(subj){
  M  = meta[["heat_crating"]][[subj]][[1]][1:36,1:186]
  Yt = meta[["ratings"]][[subj]][[1]]
  X  = meta[["heat_level"]][[subj]][[1]][1:36]
  
  colnames(M) = paste0('M',1:ncol(M))
  
  Y = rep(NA, 36)
  Y[1:length(Yt)] = Yt # Assuming the observed values represent the first trials
  
  cc = complete.cases(cbind(X,Y,M))
  
  n_miss = sum(cc)
  
  dta = data.frame(ID = subj, n_miss = n_miss, X, Y, M)[cc,]
  return(dta)
}
```


```r
big_dta = get_data(1)
for(i in 2:93) big_dta = rbind(big_dta, get_data(i)) 
#Removing outliers
big_dta = big_dta[big_dta$M186 <2,]
#Setting heat level to baseline
big_dta$X = big_dta$X - min(big_dta$X)
```


# Applying functional mediation


```r
fbootstrap_ML <- function(dta, index){
  dta = dta[index,]
  Y = dta[,"Y"]
  X = dta[,"X"]
  M = dta[,grep('M', names(dta))]
  
  result = sfs_Mediation(X,Y,t(M),
                         mediatorMethod="fosr2s",
                         splinepars_fosr2s=list(nbasis = 15, norder = 4, 
                                                basistype = "bspline"), 
                         outcomeMethod="fgam",
                         splinepars_fgam=list(bs="ps",m=c(2, 1), k=15),
                         plot=FALSE, boot=TRUE)
  return(result)
}
```

Setting seed:

```r
if(cluster) RNGkind("L'Ecuyer-CMRG")
seed = 340953
set.seed(seed)
```


```r
X = big_dta$X
Y = big_dta$Y
M = big_dta[,grep('M', colnames(big_dta))]

nbasis = 50
norder = 6
lambda = 1e-10 #1e-14 when 0.01 sd_epsilon
plots = FALSE
```


```r
time1 = Sys.time()
mediation = sfs_Mediation(x = X,y = Y,m = t(M),mediatorMethod="fosr2s",
                     splinepars_fosr2s=list(nbasis = 15, norder = 4, 
                                            basistype = "bspline"), 
                     outcomeMethod="fgam",
                     splinepars_fgam=list(bs="ps",m=c(2, 1), k=15),plot=FALSE, 
                     boot=FALSE)
time2 = Sys.time()

time3 = Sys.time()
stat = boot(data = big_dta, statistic = fbootstrap_ML,R = 510,sim = "ordinary",stype = "i",parallel = "multicore",ncpus=20)
time4 = Sys.time()
time4 - time3
```

```
## Time difference of 17.1872 mins
```

```r
save(mediation,stat, file=paste0("Data_analysis_cPain_15_",format(Sys.time(), "%Y%m%d-%H%M"),".Rdata"))
```
