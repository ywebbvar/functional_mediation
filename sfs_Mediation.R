#' Functional Mediation. Scalar stimulus - functional mediator - scalar outcome.
#' 
#' Fits a functional mediation model. Uses fda and refund packages
#'
#' Fits a functional model as in Lindquist 2012 JASA, but the defaults are different from that paper.
#' @param x a numeric vector with independent variable (treatment assignment)
#' @param y a numeric vector with final outcomes
#' @param m a (\code{T} by \code{N}) matrix with mediator values. Each columns represents one observed functional mediator.
#' @param mediatorMethod the method for model for mediator. "fRegress" fits a model by penalized from \code{fda} package, and  "fosr2" uses two-step function-on-scalar regression from \code{refund} package
#' @param outcomeMethod functional regression method to use on outcome. "fgam" defaults to quartic P-splines with 2nd derivative order penalty and at most 50-dimensional basis  
#' @param nbasis an integer variable specifying the number of basis functions. Argument for \code{\link[fda]{create.bspline.basis}}
#' @param norder an integer specifying the order of b-splines, which is one higher than their degree. Argument for \code{\link[fda]{create.bspline.basis}}
#' @param lambda a nonnegative real number specifying the amount of smoothing to be applied to the estimated functional parameters.
#' @param pen a nonnegative real number specifying the amount of penalization to be applied to the estimated functional parameter using an instrument.
#' @param plot is a logical scalar, if TRUE, the function will produce plots of the a, b and ab paths.
#' @param boot is a logical scalar, if TRUE, the function will only output the ab path, useful when bootstrap it.
#' @return afun A functional regression coefficient corresponding to a path
#' @return bfun A functional regression coefficient corresponding to b path
#' @return ab The ab effect
#' @return c The c effect
#' @return cp The c' effect 
#' @author Martin Lindquist <mlindquist@@jhu.edu>, Yenny Webb-Vargas <yennywebb@gmail.com>
#' @examples
#' sfs_Mediation(x,y,m,...)

sfs_Mediation <- function(x,y,m,mediatorMethod="fosr2s", splinepars_fosr2s=list(nbasis = 15, norder = 4, basistype = "bspline"), outcomeMethod="fgam", splinepars_fgam=list(bs="ps",m=c(3,2)), splinepars_fRegress=list(nbasis=15,norder=4,lambda=1e-8), plot=FALSE, boot=FALSE, return_fits=FALSE){
  
  require(refund)
  require(mgcv)
  
  len   = dim(m)[1]
  N     = dim(m)[2]
  T_sup = 1
  
  timevec = seq(0,T_sup, length.out=len)
  
  # Create bspline basis set
  basis = create.bspline.basis(rangeval = c(0,T_sup), nbasis = splinepars_fRegress$nbasis, splinepars_fRegress$norder)
  
  ######################
  #### Path x -> m #####
  ######################
  
  # Convert time series data to functional data
  mfd    = Data2fd(y = m,basisobj = basis,argvals = timevec)
  mfdPar = fdPar(mfd)
  
  # Create Design Matrix
  conbas = create.constant.basis(c(0,T_sup))
  confd  = fd(matrix(1,nrow=1,ncol=N), conbas)
  
  xfdcell      = list()
  xfdcell[[1]] = confd 
  xfdcell[[2]] = x
  
  # Create basis set for beta functions
  betacell      = list()
  betabasis     = create.bspline.basis(c(0,T_sup), splinepars_fRegress$nbasis, splinepars_fRegress$norder)
  betafd1       = fd(matrix(0,nrow=splinepars_fRegress$nbasis, ncol=1), betabasis)
  betacell[[1]] = fdPar(betafd1)
  betafdj       = fd(matrix(0,nrow=splinepars_fRegress$nbasis, ncol=1), betabasis)
  betacell[[2]] = fdPar(betafdj)
  
  if(mediatorMethod=="fRegress"){
    # Solve least-squares equation
    fRegressCell_M = fRegress(mfdPar, xfdcell, betacell)
    betaestcell  = fRegressCell_M[[4]]
    afun         = betaestcell[[2]]$fd
    d1fun        = betaestcell[[1]]$fd
    
    tfine = seq(0,T_sup, length.out=len)
    af    = eval.fd(tfine,afun)
    a     = sum(af)*(tfine[2] - tfine[1])
    d1f   = eval.fd(tfine,d1fun)
    
    ResM  = (m - eval.fd(tfine, fRegressCell_M[[5]]$fd))
    
    # Calculate Standard Error
    errmat     = m - eval.fd(tfine, fRegressCell_M[[5]]$fd)
    Sigma      = errmat%*%t(errmat)/N  # Originally, there was a 20 here, I assume it is the number of observations N
    DfdPar     = fdPar(basis, 0, 1)
    y2cMap     = smooth.basis(timevec,m,DfdPar)$y2cMap
    stderrCell = fRegress.stderr(fRegressCell_M, y2cMap, Sigma)
    tmp        = stderrCell[[1]]
    
    a_stderr = eval.fd(tfine, tmp[[2]]) # Std Error of a-function
    d1_stderr= eval.fd(tfine, tmp[[1]]) # Std Error of d1-function
    
    #plotbeta(betaestlist = fRegressCell_M$betaestlist, betastderrlist = stderrCell$betastderrlist)
  }else{
    if(mediatorMethod=="fosr2s"){
      tfine = seq(0,T_sup, length.out=len)
      
      fit_m = fosr2s(Y = t(m), cbind(int=1,x=x), argvals = tfine, nbasis = splinepars_fosr2s$nbasis, norder = splinepars_fosr2s$norder, basistype = splinepars_fosr2s$basistype)
      af  = fit_m$est.func[,2]
      a   = sum(af)*(tfine[2] - tfine[1])
      d1f = fit_m$est.func[,1]
      
      ResM     = t(t(m) - fit_m$yhat)
      a_stderr = fit_m$se.func[,2]
      d1_stderr= fit_m$se.func[,1]
    }else{stop("mediatorMethod must be 'fRegress' or 'fosr2s'")}
  }
  
  ######################
  ## Path x, m -> y  ###
  ######################
  
  if(outcomeMethod == "fRegress"){
    
    # Create Design matrix
    mfdcell      = list()
    
    mfdcell[[1]] = confd
    basis        = create.bspline.basis(c(0,T_sup), splinepars_fRegress$nbasis, splinepars_fRegress$norder)
    mfdcell[[2]] = Data2fd(argvals = timevec, y = m, basis)
    mfdcell[[3]] = fd(matrix(x,nrow=1,ncol=N),conbas)
    
    # Response variable
    yfdPar = y
    
    # Create basis set for beta functions
    betacell      = list()
    betafd1       = fd(1,conbas)
    betacell[[1]] = fdPar(betafd1)
    betafdj       = fd(rep(0,splinepars_fRegress$nbasis), basis)
    betafdPar     = fdPar(betafdj, lambda=splinepars_fRegress$lambda)
    betacell[[2]] = betafdPar
    betacell[[3]] = fdPar(betafd1)
    
    # Solve least-squares equation
    fRegressCell_Y = fRegress(yfdPar, mfdcell, betacell)
    betaestcell  = fRegressCell_Y[[4]]
    bfun         = betaestcell[[2]]$fd
    
    # Calculate Standard Error
    errmat     = y - fRegressCell_Y[[5]]
    ResY       = errmat
    Sigma      = errmat%*%t(errmat)/N  # Originally, there was a 20 here, I assume it is the number of observations N
    DfdPar     = fdPar(basis, 0, 1)
    y2cMap     = smooth.basis(timevec,m,DfdPar)$y2cMap
    stderrCell = fRegress.stderr(fRegressCell_Y, y2cMap, Sigma)
    tmp        = stderrCell[[1]]
    
    b_stderr  = eval.fd(tfine, tmp[[2]]) # Std Error of b-function
    g_stderr  = tmp[[3]]$coefs # Std Error of gamma
    d2_stderr = tmp[[1]]$coefs # Std Error of intercept d2
    
        
    bf = eval.fd(tfine,bfun)            # Evaluate b-function  
    b  = sum(bf)*(tfine[2]-tfine[1])    # Integral of b-function: b = \int bf(t) dt gives b-path
    
    if(mediatorMethod=="fRegress"){ 
      abf = eval.fd(tfine,afun*bfun)      # ab-functions
    }else{
      abf = af*bf
    }
    ab  = sum(abf)*(tfine[2]-tfine[1])  # Integral of ab-function: ab = \int af(t)bf(t) dt gives ab-path
    
    g  = eval.fd(tfine, betaestcell[[3]]$fd)[1] # gamma (direct effect)
    
    tmp = eval.fd(tfine, betaestcell[[1]]$fd) # d2 intercept
    d2  = tmp[1]
    
  }else{
    if(outcomeMethod=="fgam"){
      if(is.null(splinepars_fgam[["k"]])) splinepars_fgam[["k"]] = ifelse(N < 52, N-2, 50)
      m = t(m)
      fit  = fgam(y ~ x + lf(m,splinepars=splinepars_fgam)) # Defaults to quartic (m[1]=3) P-splines (bs="ps") with 2nd derivative order penalty (m[2]=2), and at most 50-dimensional basis 
      bfun = fit$smooth[[1]]
      P   = est_se_fgam(fit, term=1,n=len)
      bf  = P$estimate
      b   = sum(bf)*(tfine[2] - tfine[1])
      
      ResY     = fit$residuals
      b_stderr = P$se
      
      abf = af*bf
      ab  = sum(abf)*(tfine[2]-tfine[1])  # Integral of ab-function: ab = \int af(t)bf(t) dt gives ab-path
      
      g  = fit$coef["x"]
      d2 = fit$coef["(Intercept)"]
      
      g_stderr  = summary(fit)$se["x"]
      d2_stderr = summary(fit)$se["(Intercept)"]
    }else{stop("outcomeMethod must be 'fgam' or 'fRegress'")}
  }
  
  # Plot results
  if(plot==TRUE){
    mipar = par()$mfrow
    plot(tfine, d1f, type="l", main="'d1' function")
    lines(tfine, d1f + 2*d1_stderr, col="green")
    lines(tfine, d1f - 2*d1_stderr, col="green")
        
    par(mfrow=c(3,1))
    plot(tfine, af, type="l", main="'a' function")
    lines(tfine, af + 2*a_stderr, col="green")
    lines(tfine, af - 2*a_stderr, col="green")
    
    plot(tfine, bf, type="l", main="'b' function")
    lines(tfine, bf + 2*b_stderr, col="green")
    lines(tfine, bf - 2*b_stderr, col="green")
    
    plot(tfine, abf, type="l", main="'ab' function")
    
    par(mfrow=mipar)
  }
  #plot(x,y)
  
  if(boot==TRUE) {
    result = c(abf,ab, af, bf, g, d1f,d2)
    names(result) = c(paste0('abfunction_', 1:length(abf)),
                      "ab",
                      paste0('afunction_',  1:length(af)),
                      paste0('bfunction_',   1:length(bf)),
                      "g",
                      paste0('d1function_',   1:length(d1f)),
                      "d2")
  }else{if(return_fits){result = list()
                        result[["estimates"]] = c(abf, ab, af, bf, g, d1f,d2)
                        names(result[["estimates"]]) = c(paste0('abfunction_', 1:length(abf)),
                                                         "ab",
                                                         paste0('afunction_',  1:length(af)),
                                                         paste0('bfunction_',   1:length(c(t(bf)))),
                                                         "g",
                                                         paste0('d1function_',   1:length(d1f)),
                                                         "d2")
                        result[["output_model_fit"]] = if(outcomeMethod=="fgam"){result[["output_model_fit"]] <- fit}else{result[["output_model_fit"]] <- fRegressCell_Y}
                        result[["mediator_model_fit"]] = ifelse(mediatorMethod=="fosr2s", fit_m,fRegressCell_M)
  }else{
    result = list('afunction'  = af, 'a_stderr'  = a_stderr, 
                  'd1function' = d1f,'d1_stderr' = d1_stderr,  
                  'bfunction'  = bf, 'b_stderr'  = b_stderr, 
                  'g'          = g,  'g_stderr'  = g_stderr,
                  'd2'         = d2, 'd2_stderr' = d2_stderr,  
                  'abfunction' = abf, 
                  "ab" = ab,
                  'x' = x, 
                  'y' = y, 
                  'm' = m,
                  'tfine' = tfine,'ResM'= ResM,'ResY'= ResY)
    
  }}
  
  
  return(result)
  
}