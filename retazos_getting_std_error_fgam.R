# Finding Nemo (the standard deviation of a GAM object) y is scalar vector, x is
# scalar vector, m is a matrix with functional covariates by row. It has two 
# extra functions from plot functions. The `sub.edf` is within `plot.gam`, and 
# is not really needed. The `myplot.mgcv.smooth` is my copy of the
# `plot.mgcv.smooth` that exists somewhere but is not saved nor has a help file
# anywhere. This last one computes a matrix that is used to predict the
# functional parameter.

sub.edf <- function(lab, edf) {
  pos <- regexpr(":", lab)[1]
  if (pos < 0) {
    pos <- nchar(lab) - 1
    lab <- paste(substr(lab, start = 1, stop = pos), 
                 ",", round(edf, digits = 2), ")", sep = "")
  }
  else {
    lab1 <- substr(lab, start = 1, stop = pos - 2)
    lab2 <- substr(lab, start = pos - 1, stop = nchar(lab))
    lab <- paste(lab1, ",", round(edf, digits = 2), lab2, 
                 sep = "")
  }
  lab
}

myplot.mgcv.smooth <- function (x, P = NULL, data = NULL, label = "", se1.mult = 1, 
                                se2.mult = 2, partial.resids = FALSE, rug = TRUE, se = TRUE, 
                                scale = -1, n = 100, n2 = 40, pers = FALSE, theta = 30, phi = 30, 
                                jit = FALSE, xlab = NULL, ylab = NULL, main = NULL, ylim = NULL, 
                                xlim = NULL, too.far = 0.1, shade = FALSE, shade.col = "gray80", 
                                shift = 0, trans = I, by.resids = FALSE, scheme = 0, ...) 
{
  sp.contour <- function(x, y, z, zse, xlab = "", ylab = "", 
                         zlab = "", titleOnly = FALSE, se.plot = TRUE, se.mult = 1, 
                         trans = I, shift = 0, ...) {
    gap <- median(zse, na.rm = TRUE)
    zr <- max(trans(z + zse + shift), na.rm = TRUE) - min(trans(z - 
                                                                  zse + shift), na.rm = TRUE)
    n <- 10
    while (n > 1 && zr/n < 2.5 * gap) n <- n - 1
    zrange <- c(min(trans(z - zse + shift), na.rm = TRUE), 
                max(trans(z + zse + shift), na.rm = TRUE))
    zlev <- pretty(zrange, n)
    yrange <- range(y)
    yr <- yrange[2] - yrange[1]
    xrange <- range(x)
    xr <- xrange[2] - xrange[1]
    ypos <- yrange[2] + yr/10
    args <- as.list(substitute(list(...)))[-1]
    args$x <- substitute(x)
    args$y <- substitute(y)
    args$type = "n"
    args$xlab <- args$ylab <- ""
    args$axes <- FALSE
    do.call("plot", args)
    cs <- (yr/10)/strheight(zlab)
    if (cs > 1) 
      cs <- 1
    tl <- strwidth(zlab)
    if (tl * cs > 3 * xr/10) 
      cs <- (3 * xr/10)/tl
    args <- as.list(substitute(list(...)))[-1]
    n.args <- names(args)
    zz <- trans(z + shift)
    args$x <- substitute(x)
    args$y <- substitute(y)
    args$z <- substitute(zz)
    if (!"levels" %in% n.args) 
      args$levels <- substitute(zlev)
    if (!"lwd" %in% n.args) 
      args$lwd <- 2
    if (!"labcex" %in% n.args) 
      args$labcex <- cs * 0.65
    if (!"axes" %in% n.args) 
      args$axes <- FALSE
    if (!"add" %in% n.args) 
      args$add <- TRUE
    do.call("contour", args)
    if (is.null(args$cex.main)) 
      cm <- 1
    else cm <- args$cex.main
    if (titleOnly) 
      title(zlab, cex.main = cm)
    else {
      xpos <- xrange[1] + 3 * xr/10
      xl <- c(xpos, xpos + xr/10)
      yl <- c(ypos, ypos)
      lines(xl, yl, xpd = TRUE, lwd = args$lwd)
      text(xpos + xr/10, ypos, zlab, xpd = TRUE, pos = 4, 
           cex = cs * cm, off = 0.5 * cs * cm)
    }
    if (is.null(args$cex.axis)) 
      cma <- 1
    else cma <- args$cex.axis
    axis(1, cex.axis = cs * cma)
    axis(2, cex.axis = cs * cma)
    box()
    if (is.null(args$cex.lab)) 
      cma <- 1
    else cma <- args$cex.lab
    mtext(xlab, 1, 2.5, cex = cs * cma)
    mtext(ylab, 2, 2.5, cex = cs * cma)
    if (!"lwd" %in% n.args) 
      args$lwd <- 1
    if (!"lty" %in% n.args) 
      args$lty <- 2
    if (!"col" %in% n.args) 
      args$col <- 2
    if (!"labcex" %in% n.args) 
      args$labcex <- cs * 0.5
    zz <- trans(z + zse + shift)
    args$z <- substitute(zz)
    do.call("contour", args)
    if (!titleOnly) {
      xpos <- xrange[1]
      xl <- c(xpos, xpos + xr/10)
      lines(xl, yl, xpd = TRUE, lty = args$lty, col = args$col)
      text(xpos + xr/10, ypos, paste("-", round(se.mult), 
                                     "se", sep = ""), xpd = TRUE, pos = 4, cex = cs * 
             cm, off = 0.5 * cs * cm)
    }
    if (!"lty" %in% n.args) 
      args$lty <- 3
    if (!"col" %in% n.args) 
      args$col <- 3
    zz <- trans(z - zse + shift)
    args$z <- substitute(zz)
    do.call("contour", args)
    if (!titleOnly) {
      xpos <- xrange[2] - xr/5
      xl <- c(xpos, xpos + xr/10)
      lines(xl, yl, xpd = TRUE, lty = args$lty, col = args$col)
      text(xpos + xr/10, ypos, paste("+", round(se.mult), 
                                     "se", sep = ""), xpd = TRUE, pos = 4, cex = cs * 
             cm, off = 0.5 * cs * cm)
    }
  }
  if (is.null(P)) {
    if (!x$plot.me || x$dim > 2) 
      return(NULL)
    if (x$dim == 1) {
      raw <- data[x$term][[1]]
      if (is.null(xlim)) 
        xx <- seq(min(raw), max(raw), length = n)
      else xx <- seq(xlim[1], xlim[2], length = n)
      if (x$by != "NA") {
        by <- rep(1, n)
        dat <- data.frame(x = xx, by = by)
        names(dat) <- c(x$term, x$by)
      }
      else {
        dat <- data.frame(x = xx)
        names(dat) <- x$term
      }
      X <- PredictMat(x, dat)
      if (is.null(xlab)) 
        xlabel <- x$term
      else xlabel <- xlab
      if (is.null(ylab)) 
        ylabel <- label
      else ylabel <- ylab
      if (is.null(xlim)) 
        xlim <- range(xx)
      return(list(X = X, x = xx, scale = TRUE, se = TRUE, 
                  raw = raw, xlab = xlabel, ylab = ylabel, main = main, 
                  se.mult = se1.mult, xlim = xlim))
    }
    else {
      xterm <- x$term[1]
      if (is.null(xlab)) 
        xlabel <- xterm
      else xlabel <- xlab
      yterm <- x$term[2]
      if (is.null(ylab)) 
        ylabel <- yterm
      else ylabel <- ylab
      raw <- data.frame(x = as.numeric(data[xterm][[1]]), 
                        y = as.numeric(data[yterm][[1]]))
      n2 <- max(10, n2)
      if (is.null(xlim)) 
        xm <- seq(min(raw$x), max(raw$x), length = n2)
      else xm <- seq(xlim[1], xlim[2], length = n2)
      if (is.null(ylim)) 
        ym <- seq(min(raw$y), max(raw$y), length = n2)
      else ym <- seq(ylim[1], ylim[2], length = n2)
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
      else {
        dat <- data.frame(x = xx, y = yy)
        names(dat) <- c(xterm, yterm)
      }
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
    }
  }
  else {
    if (se) {
      if (x$dim == 1) {
        if (scheme == 1) 
          shade <- TRUE
        ul <- P$fit + P$se
        ll <- P$fit - P$se
        if (scale == 0 && is.null(ylim)) {
          ylimit <- c(min(ll), max(ul))
          if (partial.resids) {
            max.r <- max(P$p.resid, na.rm = TRUE)
            if (max.r > ylimit[2]) 
              ylimit[2] <- max.r
            min.r <- min(P$p.resid, na.rm = TRUE)
            if (min.r < ylimit[1]) 
              ylimit[1] <- min.r
          }
        }
        if (!is.null(ylim)) 
          ylimit <- ylim
        if (shade) {
          plot(P$x, trans(P$fit + shift), type = "n", 
               xlab = P$xlab, ylim = trans(ylimit + shift), 
               xlim = P$xlim, ylab = P$ylab, main = P$main, 
               ...)
          polygon(c(P$x, P$x[n:1], P$x[1]), trans(c(ul, 
                                                    ll[n:1], ul[1]) + shift), col = shade.col, 
                  border = NA)
          lines(P$x, trans(P$fit + shift), ...)
        }
        else {
          plot(P$x, trans(P$fit + shift), type = "l", 
               xlab = P$xlab, ylim = trans(ylimit + shift), 
               xlim = P$xlim, ylab = P$ylab, main = P$main, 
               ...)
          if (is.null(list(...)[["lty"]])) {
            lines(P$x, trans(ul + shift), lty = 2, ...)
            lines(P$x, trans(ll + shift), lty = 2, ...)
          }
          else {
            lines(P$x, trans(ul + shift), ...)
            lines(P$x, trans(ll + shift), ...)
          }
        }
        if (partial.resids && (by.resids || x$by == "NA")) {
          if (length(P$raw) == length(P$p.resid)) {
            if (is.null(list(...)[["pch"]])) 
              points(P$raw, trans(P$p.resid + shift), 
                     pch = ".", ...)
            else points(P$raw, trans(P$p.resid + shift), 
                        ...)
          }
          else {
            warning("Partial residuals do not have a natural x-axis location for linear functional terms")
          }
        }
        if (rug) {
          if (jit) 
            rug(jitter(as.numeric(P$raw)), ...)
          else rug(as.numeric(P$raw), ...)
        }
      }
      else if (x$dim == 2) {
        P$fit[P$exclude] <- NA
        if (pers) 
          scheme <- 1
        if (scheme == 1) {
          persp(P$x, P$y, matrix(trans(P$fit + shift), 
                                 n2, n2), xlab = P$xlab, ylab = P$ylab, zlab = P$main, 
                ylim = P$ylim, xlim = P$xlim, theta = theta, 
                phi = phi, ...)
        }
        else if (scheme == 2) {
          image(P$x, P$y, matrix(trans(P$fit + shift), 
                                 n2, n2), xlab = P$xlab, ylab = P$ylab, main = P$main, 
                xlim = P$xlim, ylim = P$ylim, col = heat.colors(50), 
                ...)
          contour(P$x, P$y, matrix(trans(P$fit + shift), 
                                   n2, n2), add = TRUE, col = 3, ...)
          if (rug) {
            if (is.null(list(...)[["pch"]])) 
              points(P$raw$x, P$raw$y, pch = ".", ...)
            else points(P$raw$x, P$raw$y, ...)
          }
        }
        else {
          sp.contour(P$x, P$y, matrix(P$fit, n2, n2), 
                     matrix(P$se, n2, n2), xlab = P$xlab, ylab = P$ylab, 
                     zlab = P$main, titleOnly = !is.null(main), 
                     se.mult = 1, trans = trans, shift = shift, 
                     ...)
          if (rug) {
            if (is.null(list(...)[["pch"]])) 
              points(P$raw$x, P$raw$y, pch = ".", ...)
            else points(P$raw$x, P$raw$y, ...)
          }
        }
      }
      else {
        warning("no automatic plotting for smooths of more than two variables")
      }
    }
    else {
      if (x$dim == 1) {
        if (scale == 0 && is.null(ylim)) {
          if (partial.resids) 
            ylimit <- range(P$p.resid, na.rm = TRUE)
          else ylimit <- range(P$fit)
        }
        if (!is.null(ylim)) 
          ylimit <- ylim
        plot(P$x, trans(P$fit + shift), type = "l", xlab = P$xlab, 
             ylab = P$ylab, ylim = trans(ylimit + shift), 
             xlim = P$xlim, main = P$main, ...)
        if (rug) {
          if (jit) 
            rug(jitter(as.numeric(P$raw)), ...)
          else rug(as.numeric(P$raw), ...)
        }
        if (partial.resids && (by.resids || x$by == "NA")) {
          if (is.null(list(...)[["pch"]])) 
            points(P$raw, trans(P$p.resid + shift), pch = ".", 
                   ...)
          else points(P$raw, trans(P$p.resid + shift), 
                      ...)
        }
      }
      else if (x$dim == 2) {
        P$fit[P$exclude] <- NA
        if (!is.null(main)) 
          P$title <- main
        if (pers) 
          scheme <- 1
        if (scheme == 1) {
          persp(P$x, P$y, matrix(trans(P$fit + shift), 
                                 n2, n2), xlab = P$xlab, ylab = P$ylab, zlab = P$main, 
                theta = theta, phi = phi, xlim = P$xlim, 
                ylim = P$ylim, ...)
        }
        else if (scheme == 2) {
          image(P$x, P$y, matrix(trans(P$fit + shift), 
                                 n2, n2), xlab = P$xlab, ylab = P$ylab, main = P$main, 
                xlim = P$xlim, ylim = P$ylim, col = heat.colors(50), 
                ...)
          contour(P$x, P$y, matrix(trans(P$fit + shift), 
                                   n2, n2), add = TRUE, col = 3, ...)
          if (rug) {
            if (is.null(list(...)[["pch"]])) 
              points(P$raw$x, P$raw$y, pch = ".", ...)
            else points(P$raw$x, P$raw$y, ...)
          }
        }
        else {
          contour(P$x, P$y, matrix(trans(P$fit + shift), 
                                   n2, n2), xlab = P$xlab, ylab = P$ylab, main = P$main, 
                  xlim = P$xlim, ylim = P$ylim, ...)
          if (rug) {
            if (is.null(list(...)[["pch"]])) 
              points(P$raw$x, P$raw$y, pch = ".", ...)
            else points(P$raw$x, P$raw$y, ...)
          }
        }
      }
      else {
        warning("no automatic plotting for smooths of more than one variable")
      }
    }
  }
}
################## 
# My code
##################

fit = fgam(y ~ x + lf(m,splinepars=list(bs="ps", k=ifelse(N < 52, N-2, 50),m=c(3,2)))) # Defaults to quartic (m[1]=3) P-splines (bs="ps") with 2nd derivative order penalty (m[2]=2), and at most 50-dimensional basis 

predictions = predict(fit)
newdata = data.frame(m.tmat = seq(0,1,length=len), L.m = seq(1,1,length=len))
sm = fit$smooth[[1]]
bf = PredictMat(sm, newdata)%*%fit$coef[-(1:2)]
b   = sum(bf)*(tfine[2] - tfine[1])

ResY     = fit$residuals


smooth_term = fit$smooth[[1]]

first <- smooth_term$first.para
last <- smooth_term$last.para

edf <- sum(fit$edf[first:last])
term.lab <- sub.edf(smooth_term$label, edf)
attr(smooth_term, "coefficients") <- fit$coefficients[first:last]

raw <- fit$model[smooth_term$term][[1]]
n = 100
xx <- seq(min(raw), max(raw), length = n)
by  <- rep(1, n)
dat <- data.frame(x = seq(min(raw), max(raw), length = n), by = rep(1, n))
names(dat) <- c(smooth_term$term, smooth_term$by)
P = list()
P$X <- PredictMat(smooth_term, dat)
P$se.mult = qnorm(0.975)


P = myplot.mgcv.smooth(fit$smooth[[1]], P = NULL, data = fit$model)#, 

se.fit <- sqrt(pmax(0, rowSums((P$X %*% 
                                  fit$Vp[first:last, first:last, drop = FALSE]) * 
                                 P$X)))
P$se <- se.fit * P$se.mult

b_stderr = P$se

