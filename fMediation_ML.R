#' Functional Mediation
#' 
#' Fits a functional mediation model. Uses fda package
#'
#' Fits a functional model as in Lindquist 2012 JASA.
#' @param x A 1 by N vector with treatment assignment
#' @param y A 1 by N vector with outcome
#' @param m A len(time) by N (repetitions) vector with mediator values
#' @param nbasis
#' @param norder
#' @param estimate
#' @param lambda
#' @param pen 
#' @return afun A functional regression coefficient corresponding to a path
#' @return bfun A functional regression coefficient corresponding to b path
#' @return ab The ab effect
#' @return c The c effect
#' @return cp The c' effect 
#' @author Martin Lindquist <mlindquist@@jhu.edu>
#' @examples
#' fMediation_ML(x,y,m,...)

fMediation_ML <- function(x,y,m,nbasis,norder,estimate=2,lambda=10,pen=0.1){
  
}