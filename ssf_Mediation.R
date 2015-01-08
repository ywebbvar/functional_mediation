#' Functional Mediation. Scalar stimulus - scalar mediator - functional outcome.
#' 
#' Fits a scalar-scalar-functional mediation model. Uses fda and refund packages
#'
#' @param x a numeric vector with independent variable (treatment assignment) values
#' @param m a numeric vector with mediator values
#' @param y a (\code{T} by \code{N}) matrix with final outcome values. Each column represents one observed functional outcome.
#' @param outcomeMethod the method for model for mediator. "fRegress" fits a model by penalized from \code{fda} package, and  "fosr2s" uses two-step function-on-scalar regression from \code{refund} package
#' @param nbasis an integer variable specifying the number of basis functions. Argument for \code{\link[fda]{create.bspline.basis}}
#' @param norder an integer specifying the order of b-splines, which is one higher than their degree. Argument for \code{\link[fda]{create.bspline.basis}}
#' @param plot is a logical scalar, if TRUE, the function will produce plots of the a, b and ab paths.
#' @param boot is a logical scalar, if TRUE, the function will only output the ab path, useful when bootstrap it.
#' @return gfun A functional regression coefficient corresponding to the gamma path
#' @return bfun A functional regression coefficient corresponding to the beta path
#' @return ab The alpha-beta effect
#' @return cfun A functional regression coefficient for the total effect at a particular timepoint
#' @author Yenny Webb-Vargas <ywebbva1@@jhu.edu> and Martin Lindquist <mlindquist@@jhu.edu>
#' @examples
#' ssf_Mediation(x,y,m,...)

ssf_Mediation <- function(x,y,m,outcomeMethod="fosr2s", nbasis,norder,plot=FALSE, boot=FALSE){
  
  require(refund)
  
  len   = dim(y)[1]
  N     = dim(y)[2]
  T_sup = 1
  
  timevec = seq(0,T_sup, length.out=len)
  
  # Create bspline basis set
  basis = create.bspline.basis(rangeval = c(0,T_sup), nbasis = nbasis, norder)
  
  ######################
  #### Path x -> m #####
  ######################
  m_model = lm(m~x)
  d1  = m_model$coef["(Intercept)"]
  a   = m_model$coef["x"]
  
  d1_stderr = summary(m_model)$coef["(Intercept)","Std. Error"]
  a_stderr  = summary(m_model)$coef["x","Std. Error"]
  
  ResM = m_model$res
  
  ######################
  ## Path x, m -> y  ###
  ######################
  
  # Convert time series data to functional data
  yfd    = Data2fd(y = y,basisobj = basis,argvals = timevec)
  yfdPar = fdPar(yfd)
  
  # Create Design Matrix
  conbas = create.constant.basis(c(0,T_sup))
  confd  = fd(matrix(1,nrow=1,ncol=N), conbas)
  
  xfdcell      = list()
  xfdcell[[1]] = confd 
  xfdcell[[2]] = x
  xfdcell[[3]] = m
  
  # Create basis set for beta functions
  betacell      = list()
  betabasis     = create.bspline.basis(c(0,T_sup), nbasis, norder)
  betafd1       = fd(matrix(0,nrow=nbasis, ncol=1), betabasis)
  betacell[[1]] = fdPar(betafd1)
  betacell[[2]] = fdPar(betafd1)
  betacell[[3]] = fdPar(betafd1)
  
  if(outcomeMethod=="fRegress"){
    # Solve least-squares equation
    fRegressCell = fRegress(yfdPar, xfdcell, betacell)
    betaestcell  = fRegressCell[[4]]
    d2fun        = betaestcell[[1]]$fd
    gfun         = betaestcell[[2]]$fd
    bfun         = betaestcell[[3]]$fd
    
    tfine = seq(0,T_sup, length.out=len)
    d2f   = eval.fd(tfine,d2fun)
    gf    = eval.fd(tfine,gfun)
    bf    = eval.fd(tfine,bfun)
    
    ResY  = (y - eval.fd(tfine, fRegressCell$yhatfdobj$fd))
    
    # Calculate Standard Error
    errmat     = y - eval.fd(tfine, fRegressCell$yhatfdobj$fd)
    Sigma      = errmat%*%t(errmat)/N 
    DfdPar     = fdPar(basis, 0, 1)
    y2cMap     = smooth.basis(timevec,y,DfdPar)$y2cMap
    stderrCell = fRegress.stderr(fRegressCell, y2cMap, Sigma)
    tmp        = stderrCell[[1]]
    
    d2_stderr = eval.fd(tfine, stderrCell$betastderrlist[[1]]) # Std Error of d2-function
    g_stderr  = eval.fd(tfine, stderrCell$betastderrlist[[2]]) # Std Error of g-function
    b_stderr  = eval.fd(tfine, stderrCell$betastderrlist[[3]]) # Std Error of b-function
    plotbeta(betaestlist = fRegressCell$betaestlist, betastderrlist = stderrCell$betastderrlist)
  }else{
    if(outcomeMethod=="fosr2s"){
      tfine = seq(0,T_sup, length.out=len)
      
      fit = fosr2s(Y = t(y), cbind(int=1,x=x,m=m), argvals = tfine, nbasis = 15, norder = 4, basistype = "bspline")
      d2f = fit$est.func[,1]
      gf  = fit$est.func[,2]
      bf  = fit$est.func[,3]
      
      ResY     = t(t(y) - fit$yhat)
      d2_stderr = fit$se.func[,1] # Std Error of d2-function
      g_stderr  = fit$se.func[,2] # Std Error of g-function
      b_stderr  = fit$se.func[,3] # Std Error of b-function
    }else{stop("outcomeMethod must be 'fRegress' or 'fosr2s'")}
  }
  
  
   abf = a*bf
  
  # Plot results
  if(plot==TRUE){
    mipar = par()$mfrow
    par(mfrow=c(2,2))
  
    plot(tfine, d2f, type="l", main="'delta_2' function")
    lines(tfine, d2f + 2*d2_stderr, col="green")
    lines(tfine, d2f - 2*d2_stderr, col="green")
  
    plot(tfine, gf, type="l", main="'gamma' function")
    lines(tfine, gf + 2*g_stderr, col="green")
    lines(tfine, gf - 2*g_stderr, col="green")
        
    plot(tfine, bf, type="l", main="'beta' function")
    lines(tfine, bf + 2*b_stderr, col="green")
    lines(tfine, bf - 2*b_stderr, col="green")
    
    plot(tfine, abf, type="l", main="'alpha*beta' function")
    par(mfrow=mipar)
  }
  #plot(x,y)
  
  m_path = list(delta1 = d1, alpha = a, delta1_stderr = d1_stderr, alpha_stderr= a_stderr, model= m_model)
  
  y_path = list(delta2_function = d2f, gamma_function = gf, beta_function=bf, delta2_stderr=d2_stderr, gamma_stderr = g_stderr, beta_stderr=b_stderr)
  
  if(boot==FALSE) {result = list(m_path=m_path, y_path = y_path, data = list('x' = x, 'y' = y, 'm' = m),tfine = tfine, ResM= ResM,ResY= ResY)
  }else{
    result = abf
  }
  
  return(result)
  
}