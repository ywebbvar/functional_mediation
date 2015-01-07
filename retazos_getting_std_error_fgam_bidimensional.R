
n2 = 40; too.far = 0.1

first <- x$smooth[[i]]$first.para
last <- x$smooth[[i]]$last.para

edf <- sum(x$edf[first:last])
term.lab <- sub.edf(x$smooth[[i]]$label, edf)
attr(x$smooth[[i]], "coefficients") <- x$coefficients[first:last]


x = fit
data = fit$model
if (!x$plot.me || x$dim > 2) 
  
  xterm <- x$term[1]
#if (is.null(xlab)) 
  xlabel <- xterm
#else xlabel <- xlab
yterm <- x$term[2]
#if (is.null(ylab)) 
  ylabel <- yterm
#else ylabel <- ylab
raw <- data.frame(x = as.numeric(data[xterm][[1]]), 
                  y = as.numeric(data[yterm][[1]]))
n2 <- max(10, n2)
#if (is.null(xlim)) 
  xm <- seq(min(raw$x), max(raw$x), length = n2)
#else xm <- seq(xlim[1], xlim[2], length = n2)
#if (is.null(ylim)) 
  ym <- seq(min(raw$y), max(raw$y), length = n2)
#else ym <- seq(ylim[1], ylim[2], length = n2)
xx <- rep(xm, n2)
yy <- rep(ym, rep(n2, n2))
if (too.far > 0) 
  exclude <- exclude.too.far(xx, yy, raw$x, raw$y, 
                             dist = too.far)
else exclude <- rep(FALSE, n2 * n2)
if (x$by != "NA") {
  by <- rep(1, n2^2)
  dat <- data.frame(x = xx, y = yy, by = by)
  names(dat) <- c(xterm, yterm, x$by)
}
#else {
#  dat <- data.frame(x = xx, y = yy)
#  names(dat) <- c(xterm, yterm)
#}
X <- PredictMat(x, dat)
if (is.null(main)) {
  main <- label
}
if (is.null(ylim)) 
  ylim <- range(ym)
if (is.null(xlim)) 
  xlim <- range(xm)
return(list(X = X, x = xm, y = ym, scale = FALSE, 
            se = TRUE, raw = raw, xlab = xlabel, ylab = ylabel, 
            main = main, se.mult = se2.mult, ylim = ylim, 
            xlim = xlim, exclude = exclude))

p <- x$coefficients[first:last]
offset <- attr(P$X, "offset")
#if (is.null(offset)) 
  P$fit <- P$X %*% p
#else P$fit <- P$X %*% p + offset
if (!is.null(P$exclude)) 
  P$fit[P$exclude] <- NA
if (se && P$se) {
#  if (seWithMean && attr(x$smooth[[i]], "nCons") > 
#        0) {
#    if (length(x$cmX) < ncol(x$Vp)) 
#      x$cmX <- c(x$cmX, rep(0, ncol(x$Vp) - length(x$cmX)))
#    X1 <- matrix(x$cmX, nrow(P$X), ncol(x$Vp), 
#                 byrow = TRUE)
#    meanL1 <- x$smooth[[i]]$meanL1
#    if (!is.null(meanL1)) 
#      X1 <- X1/meanL1
#    X1[, first:last] <- P$X
#    se.fit <- sqrt(pmax(0, rowSums((X1 %*% x$Vp) * 
#                                     X1)))
#  }
#  else{ 
    se.fit <- sqrt(pmax(0, rowSums((P$X %*% 
                                         x$Vp[first:last, first:last, drop = FALSE]) * 
                                        P$X)))
#  }
  if (!is.null(P$exclude)) 
    P$se.fit[P$exclude] <- NA
}
if (partial.resids) {
  P$p.resid <- fv.terms[, length(order) + i] + 
    w.resid
}
if (se && P$se) 
  P$se <- se.fit * P$se.mult
P$X <- NULL
P$plot.me <- TRUE
pd[[i]] <- P
rm(P)





if (m > 0) 
  for (i in 1:m) if (pd[[i]]$plot.me && (is.null(select) || 
                                           i == select)) {
    plot(x$smooth[[i]], P = pd[[i]], partial.resids = partial.resids, 
         rug = rug, se = se, scale = scale, n = n, n2 = n2, 
         pers = pers, theta = theta, phi = phi, jit = jit, 
         xlab = xlab, ylab = ylab, main = main, ylim = ylim, 
         xlim = xlim, too.far = too.far, shade = shade, 
         shade.col = shade.col, shift = shift, trans = trans, 
         by.resids = by.resids, scheme = scheme[i], ...)
  }


shift=0

debug: if (x$dim == 2) {
  P$fit[P$exclude] <- NA
  if (pers) 
    scheme <- 1
  if (scheme == 1) {
    persp(P$x, P$y, matrix(trans(P$fit + shift), n2, n2), 
          xlab = P$xlab, ylab = P$ylab, zlab = P$main, ylim = P$ylim, 
          xlim = P$xlim, theta = theta, phi = phi, ...)
  }
  else if (scheme == 2) {
    image(P$x, P$y, matrix(
      (P$fit + shift), n2, n2), 
          xlab = P$xlab, ylab = P$ylab, main = P$main, xlim = P$xlim, 
          ylim = P$ylim, col = heat.colors(50), ...)
    contour(P$x, P$y, matrix(trans(P$fit + shift), n2, n2), 
            add = TRUE, col = 3, ...)
    if (rug) {
      if (is.null(list(...)[["pch"]])) 
        points(P$raw$x, P$raw$y, pch = ".", ...)
      else points(P$raw$x, P$raw$y, ...)
    }
  }
  
  else {
  ### This is the one it uses by default)  
    sp.contour(P$x, P$y, matrix(P$fit, n2, n2), matrix(P$se, 
                                                       n2, n2), xlab = P$xlab, ylab = P$ylab, zlab = P$main, 
               titleOnly = !is.null(main), se.mult = 1, trans = trans, 
               shift = shift, ...)
    if (rug) {
      if (is.null(list(...)[["pch"]])) 
        points(P$raw$x, P$raw$y, pch = ".", ...)
      else points(P$raw$x, P$raw$y, ...)
    }
  }
} else {
  warning("no automatic plotting for smooths of more than two variables")
}