#' Finding the estimate and standard deviation of a functional parameter in GAM
#' object
#' 
#' Takes an GAM object, selects one smooth term, and computes estimated
#' functional parameter and its standard error.
#' 
#' Takes an GAM object, selects one smooth term, and computes estimated
#' functional parameter and its standard error. It takes code from
#' plot.mgcv.smooth and plot.gam
#' @param fit an GAM object
#' @param term an integer signaling the smooth term which parameters will be computed
#' @return P A list including the predicted matrix of the smooth term, the estimate, and the standard error.
#' @author Yenny Webb-Vargas <yennywebb@gmail.com>
#' @examples
#' se_fgam(fit)
est_se_fgam = function(fit, term=1, ff=FALSE, n=100){
  require(mgcv)
  smooth_term = fit$smooth[[term]]
  
  first = smooth_term$first.para
  last  = smooth_term$last.para
  
  raw = fit$model[smooth_term$term][[1]]
  if(ff == FALSE){
    dat = data.frame(x = seq(min(raw), max(raw), length = n), by = rep(1, n))
    names(dat) = c(smooth_term$term, smooth_term$by)
    P = list()
    P$X = PredictMat(smooth_term, dat)
    P$estimate = P$X%*%fit$coefficients[first:last]
    P$se = sqrt(pmax(0, rowSums((P$X %*%fit$Vp[first:last, first:last, drop = FALSE])*P$X)))
  }else{
    dat = data.frame(x = seq(min(raw), max(raw), length = n),y = seq(min(raw), max(raw), length = n), by = rep(1, n))
    names(dat) = c(smooth_term$term, smooth_term$by)
    P = list()
    P$X = PredictMat(smooth_term, dat)
    P$estimate = P$X%*%fit$coefficients[first:last]
    P$se = sqrt(pmax(0, rowSums((P$X %*%fit$Vp[first:last, first:last, drop = FALSE])*P$X)))
  }
  return(P)
}
