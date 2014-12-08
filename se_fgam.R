#' Finding the standard deviation of an GAM object
#' 
#' Takes an GAM object and computes standard errors at each time point.
#'
#' Takes an GAM object and computes standard errors at each time point. It takes code from plot.mgcv.smooth and plot.gam
#' @param fit an GAM object
#' @return P A functional regression coefficient corresponding to a path
#' @author Yenny Webb-Vargas <yennywebb@gmail.com>
#' @examples
#' se_fgam(fit)
se_fgam = function(fit){
  n = 100
  smooth_term = fit$smooth[[1]]
  
  first = smooth_term$first.para
  last  = smooth_term$last.para
  
  raw = fit$model[smooth_term$term][[1]]
  xx = seq(min(raw), max(raw), length = n)
  by  = rep(1, n)
  dat = data.frame(x = seq(min(raw), max(raw), length = n), by = rep(1, n))
  names(dat) = c(smooth_term$term, smooth_term$by)
  P = list()
  P$X = PredictMat(smooth_term, dat)
  P$se.mult = qnorm(0.975)
  P$se.fit = sqrt(pmax(0, rowSums((P$X %*% 
                                    fit$Vp[first:last, first:last, drop = FALSE]) * 
                                   P$X)))
  P$CI = P$se.fit * P$se.mult
  return(P)
}
