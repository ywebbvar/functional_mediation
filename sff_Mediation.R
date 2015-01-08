#' Functional Mediation. Scalar stimulus - functional mediator - functional outcome.
#' 
#' Fits a functional mediation model. Uses fda and refund packages
#'
#' Fits a functional model as in Lindquist 2012 JASA.
#' @param x a numeric vector with independent variable (treatment assignment)
#' @param y a numeric vector with final outcomes
#' @param m a (\code{T} by \code{N}) matrix with mediator values. Each columns represents one observed functional mediator.
#' @param mediatorMethod the method for model for mediator. "fRegress" fits a model by penalized from \code{fda} package, and  "fosr2" uses two-step function-on-scalar regression from \code{refund} package
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
#' @author Martin Lindquist <mlindquist@@jhu.edu>
#' @examples
#' fMediation_ML(x,y,m,...)

sff_Mediation <- function(x,y,m,mediatorMethod="fosr2s", nbasis,norder,lambda=1e-8,pen=0.1, plot=FALSE, boot=FALSE){
  
  require(refund)
  
  len   = dim(m)[1]
  N     = dim(m)[2]
  T_sup = 1
  
  timevec = seq(0,T_sup, length.out=len)
  
  # Create bspline basis set
  basis = create.bspline.basis(rangeval = c(0,T_sup), nbasis = nbasis, norder)
  
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
  betabasis     = create.bspline.basis(c(0,T_sup), nbasis, norder)
  betafd1       = fd(matrix(0,nrow=nbasis, ncol=1), betabasis)
  betacell[[1]] = fdPar(betafd1)
  betafdj       = fd(matrix(0,nrow=nbasis, ncol=1), betabasis)
  betacell[[2]] = fdPar(betafdj)
  
  if(mediatorMethod=="fRegress"){
    # Solve least-squares equation
    fRegressCell = fRegress(mfdPar, xfdcell, betacell)
    betaestcell  = fRegressCell[[4]]
    afun         = betaestcell[[2]]$fd
    d1fun        = betaestcell[[1]]$fd
    
    tfine = seq(0,T_sup, length.out=len)
    af    = eval.fd(tfine,afun)
    a     = sum(af)*(tfine[2] - tfine[1])
    d1f   = eval.fd(tfine,d1fun)
    
    ResM  = (m - eval.fd(tfine, fRegressCell[[5]]$fd))
    
    # Calculate Standard Error
    errmat     = m - eval.fd(tfine, fRegressCell[[5]]$fd)
    Sigma      = errmat%*%t(errmat)/N  # Originally, there was a 20 here, I assume it is the number of observations N
    DfdPar     = fdPar(basis, 0, 1)
    y2cMap     = smooth.basis(timevec,m,DfdPar)$y2cMap
    stderrCell = fRegress.stderr(fRegressCell, y2cMap, Sigma)
    tmp        = stderrCell[[1]]
    
    a_stderr = eval.fd(tfine, tmp[[2]]) # Std Error of a-function
    d1_stderr= eval.fd(tfine, tmp[[1]]) # Std Error of d1-function
    
    #plotbeta(betaestlist = fRegressCell$betaestlist, betastderrlist = stderrCell$betastderrlist)
  }else{
    if(mediatorMethod=="fosr2s"){
      tfine = seq(0,T_sup, length.out=len)
      
      fit = fosr2s(Y = t(m), cbind(int=1,x=x), argvals = tfine, nbasis = 15, norder = 4, basistype = "bspline")
      af  = fit$est.func[,2]
      a   = sum(af)*(tfine[2] - tfine[1])
      d1f = fit$est.func[,1]
      
      ResM     = t(t(m) - fit$yhat)
      a_stderr = fit$se.func[,2]
      d1_stderr= fit$se.func[,1]
    }else{stop("mediatorMethod must be 'fRegress' or 'fosr2s'")}
  }
  
  ######################
  ## Path x, m -> y  ###
  ######################
  tfine = seq(0,T_sup, length.out=len)
    
  m = t(m)
  fit  = pffr(y ~ x + ff(m,limits="s<t", integration="riemann", splinepars=list(bs="pss", k=ifelse(N < 52, N-2, 50),m=c(3,2)))) # Defaults to quartic (m[1]=3) P-splines (bs="ps") with 2nd derivative order penalty (m[2]=2), and at most 50-dimensional basis 
    
  d2fun = fit$smooth[["s(yindex.vec)"]]
  Pd2   = est_se_fgam(fit, term="s(yindex.vec)",n=len)
  d2f   = Pd2$estimate
  d2_stderr = Pd2$se  
  
  gfun = fit$smooth[["s(yindex.vec):x"]]
  Pg   = est_se_fgam(fit, term="s(yindex.vec):x",n=len)
  gf   = Pg$estimate
  g_stderr = Pg$se  
  
  bfun = fit$smooth[["te(m.smat,m.tmat):L.m"]]
  Pb   = est_se_fgam(fit, term="te(m.smat,m.tmat):L.m",n=len, ff=TRUE, n2=len)
  bf   = Pb$estimate
  b   = sum(bf)*(tfine[2] - tfine[1])
  b_stderr = Pb$se  
  
  ResY     = fit$residuals
  
  abf = matrix(af, byrow=T, ncol=len, nrow=len)*bf
  ab  = sum(abf)*(tfine[2]-tfine[1])  # Integral of ab-function: ab(t) = \int af(t)bf(t,s)ds dt gives ab-path

  # Plot results
  if(plot==TRUE){
    mipar = par()$mfrow
    par(mfrow=c(2,1),ask = TRUE)

    plot(tfine, d1f, type="l", main="'Delta1' function")
    lines(tfine, d1f + 2*d1_stderr, col="green")
    lines(tfine, d1f - 2*d1_stderr, col="green")
    
    plot(tfine, af, type="l", main="'Alpha' function")
    lines(tfine, af + 2*a_stderr, col="green")
    lines(tfine, af - 2*a_stderr, col="green")
    
    plot(tfine, d2f, type="l", main="'Delta2' function")
    lines(tfine, d2f + 2*d2_stderr, col="green")
    lines(tfine, d2f - 2*d2_stderr, col="green")
    
    plot(tfine, gf, type="l", main="'Gamma' function")
    lines(tfine, gf + 2*g_stderr, col="green")
    lines(tfine, gf - 2*g_stderr, col="green")
  
    par(mfrow=c(1,1))
    image(t(bf), col  = gray((0:32)/32), main="Beta function")
    image(t(abf), col  = gray((0:32)/32), main="Alpha*Beta function")
    par(mfrow=mipar, ask=FALSE)
  }
    
  if(boot==FALSE) {result = list('afunction' = af, 'a' = a, 'bfunction' = bf, 'abfunction' = abf, 'b' = b, 'ab' = ab, 'cp' = cp, 'Y_intercept' = int, 'x' = x, 'y' = y, 'm' = m,'tfine' = tfine,'b_stderr' = b_stderr,'ResM'= ResM,'ResY'= ResY)
  }else{
    result = abf
  }
  
  return(result)
  
}