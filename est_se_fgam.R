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
est_se_fgam = function(fit, term=1, ff=FALSE, n=100,n2 = 40,too.far = 0.1){
  require(mgcv)
  smooth_term = fit$smooth[[term]]
  
  first = smooth_term$first.para
  last  = smooth_term$last.para
  
  raw = fit$model[smooth_term$term][[1]]
  
  x = fit
  data = fit$model
  
  if(ff == FALSE){
    dat = data.frame(x = seq(min(raw), max(raw), length = n), by = rep(1, n))
    names(dat) = c(smooth_term$term, smooth_term$by)
    P = list()
    P$X = PredictMat(smooth_term, dat)
    P$estimate = P$X%*%fit$coefficients[first:last]
    P$se = sqrt(pmax(0, rowSums((P$X %*%fit$Vp[first:last, first:last, drop = FALSE])*P$X)))
  }else{
    xterm <- fit$term[3]
    xlabel <- xterm
    yterm <- fit$term[4]
    ylabel <- yterm
    raw <- data.frame(x = as.numeric(data[xterm][[1]]), 
                      y = as.numeric(data[yterm][[1]]))
    n2 <- max(10, n2)
    xm <- seq(min(raw$x), max(raw$x), length = n2)
    ym <- seq(min(raw$y), max(raw$y), length = n2)
    xx <- rep(xm, n2)
    yy <- rep(ym, rep(n2, n2))
    if (too.far > 0) 
      exclude <- exclude.too.far(xx, yy, raw$x, raw$y, 
                                 dist = too.far)
    else exclude <- rep(FALSE, n2 * n2)
    by <- rep(1, n2^2)
    dat <- data.frame(x = xx, y = yy, by = by)
    names(dat) <- c(xterm, yterm, x$by)
    
    P = list()
    P$X = PredictMat(smooth_term, dat)
    P$estimate = P$X%*%fit$coefficients[first:last]
    P$se = sqrt(pmax(0, rowSums((P$X %*%fit$Vp[first:last, first:last, drop = FALSE])*P$X)))
  }
  return(P)
}
