> plot.gam
function (x, residuals = FALSE, rug = TRUE, se = TRUE, pages = 0, 
          select = NULL, scale = -1, n = 100, n2 = 40, pers = FALSE, 
          theta = 30, phi = 30, jit = FALSE, xlab = NULL, ylab = NULL, 
          main = NULL, ylim = NULL, xlim = NULL, too.far = 0.1, all.terms = FALSE, 
          shade = FALSE, shade.col = "gray80", shift = 0, trans = I, 
          seWithMean = FALSE, unconditional = FALSE, by.resids = FALSE, 
          scheme = 0, ...) 
{
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
  if (unconditional) {
    if (is.null(x$Vc)) 
      warning("Smoothness uncertainty corrected covariance not available")
    else x$Vp <- x$Vc
  }
  w.resid <- NULL
  if (length(residuals) > 1) {
    if (length(residuals) == length(x$residuals)) 
      w.resid <- residuals
    else warning("residuals argument to plot.gam is wrong length: ignored")
    partial.resids <- TRUE
  }
  else partial.resids <- residuals
  m <- length(x$smooth)
  if (length(scheme) == 1) 
    scheme <- rep(scheme, m)
  if (length(scheme) != m) {
    warn <- paste("scheme should be a single number, or a vector with", 
                  m, "elements")
    warning(warn)
    scheme <- rep(scheme[1], m)
  }
  order <- if (is.list(x$pterms)) 
    unlist(lapply(x$pterms, attr, "order"))
  else attr(x$pterms, "order")
  if (all.terms) 
    n.para <- sum(order == 1)
  else n.para <- 0
  if (se) {
    if (is.numeric(se)) 
      se2.mult <- se1.mult <- se
    else {
      se1.mult <- 2
      se2.mult <- 1
    }
    if (se1.mult < 0) 
      se1.mult <- 0
    if (se2.mult < 0) 
      se2.mult <- 0
  }
  else se1.mult <- se2.mult <- 1
  if (se && x$Vp[1, 1] < 0) {
    se <- FALSE
    warning("No variance estimates available")
  }
  if (partial.resids) {
    if (is.null(w.resid)) {
      if (is.null(x$residuals) || is.null(x$weights)) 
        partial.resids <- FALSE
      else {
        wr <- sqrt(x$weights)
        w.resid <- x$residuals * wr/mean(wr)
      }
    }
    if (partial.resids) 
      fv.terms <- predict(x, type = "terms")
  }
  pd <- list()
  i <- 1
  if (m > 0) 
    for (i in 1:m) {
      first <- x$smooth[[i]]$first.para
      last <- x$smooth[[i]]$last.para
      edf <- sum(x$edf[first:last])
      term.lab <- sub.edf(x$smooth[[i]]$label, edf)
      attr(x$smooth[[i]], "coefficients") <- x$coefficients[first:last]
      P <- plot(x$smooth[[i]], P = NULL, data = x$model, 
                partial.resids = partial.resids, rug = rug, se = se, 
                scale = scale, n = n, n2 = n2, pers = pers, theta = theta, 
                phi = phi, jit = jit, xlab = xlab, ylab = ylab, 
                main = main, label = term.lab, ylim = ylim, xlim = xlim, 
                too.far = too.far, shade = shade, shade.col = shade.col, 
                se1.mult = se1.mult, se2.mult = se2.mult, shift = shift, 
                trans = trans, by.resids = by.resids, scheme = scheme[i], 
                ...)
      if (is.null(P)) 
        pd[[i]] <- list(plot.me = FALSE)
      else if (is.null(P$fit)) {
        p <- x$coefficients[first:last]
        offset <- attr(P$X, "offset")
        if (is.null(offset)) 
          P$fit <- P$X %*% p
        else P$fit <- P$X %*% p + offset
        if (!is.null(P$exclude)) 
          P$fit[P$exclude] <- NA
        if (se && P$se) {
          if (seWithMean && attr(x$smooth[[i]], "nCons") > 
                0) {
            if (length(x$cmX) < ncol(x$Vp)) 
              x$cmX <- c(x$cmX, rep(0, ncol(x$Vp) - length(x$cmX)))
            X1 <- matrix(x$cmX, nrow(P$X), ncol(x$Vp), 
                         byrow = TRUE)
            meanL1 <- x$smooth[[i]]$meanL1
            if (!is.null(meanL1)) 
              X1 <- X1/meanL1
            X1[, first:last] <- P$X
            se.fit <- sqrt(pmax(0, rowSums((X1 %*% x$Vp) * 
                                             X1)))
          }
          else se.fit <- sqrt(pmax(0, rowSums((P$X %*% 
                                                 x$Vp[first:last, first:last, drop = FALSE]) * 
                                                P$X)))
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
      }
      else {
        if (partial.resids) {
          P$p.resid <- fv.terms[, length(order) + i] + 
            w.resid
        }
        P$plot.me <- TRUE
        pd[[i]] <- P
        rm(P)
      }
    }
  n.plots <- n.para
  if (m > 0) 
    for (i in 1:m) n.plots <- n.plots + as.numeric(pd[[i]]$plot.me)
  if (n.plots == 0) 
    stop("No terms to plot - nothing for plot.gam() to do.")
  if (pages > n.plots) 
    pages <- n.plots
  if (pages < 0) 
    pages <- 0
  if (pages != 0) {
    ppp <- n.plots%/%pages
    if (n.plots%%pages != 0) {
      ppp <- ppp + 1
      while (ppp * (pages - 1) >= n.plots) pages <- pages - 
        1
    }
    c <- r <- trunc(sqrt(ppp))
    if (c < 1) 
      r <- c <- 1
    if (c * r < ppp) 
      c <- c + 1
    if (c * r < ppp) 
      r <- r + 1
    oldpar <- par(mfrow = c(r, c))
  }
  else {
    ppp <- 1
    oldpar <- par()
  }
  if ((pages == 0 && prod(par("mfcol")) < n.plots && dev.interactive()) || 
        pages > 1 && dev.interactive()) 
    ask <- TRUE
  else ask <- FALSE
  if (!is.null(select)) {
    ask <- FALSE
  }
  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }
  if (scale == -1 && is.null(ylim)) {
    k <- 0
    if (m > 0) 
      for (i in 1:m) if (pd[[i]]$plot.me && pd[[i]]$scale) {
        if (se && length(pd[[i]]$se) > 1) {
          ul <- pd[[i]]$fit + pd[[i]]$se
          ll <- pd[[i]]$fit - pd[[i]]$se
          if (k == 0) {
            ylim <- c(min(ll, na.rm = TRUE), max(ul, 
                                                 na.rm = TRUE))
            k <- 1
          }
          else {
            if (min(ll, na.rm = TRUE) < ylim[1]) 
              ylim[1] <- min(ll, na.rm = TRUE)
            if (max(ul, na.rm = TRUE) > ylim[2]) 
              ylim[2] <- max(ul, na.rm = TRUE)
          }
        }
        else {
          if (k == 0) {
            ylim <- range(pd[[i]]$fit, na.rm = TRUE)
            k <- 1
          }
          else {
            if (min(pd[[i]]$fit, na.rm = TRUE) < ylim[1]) 
              ylim[1] <- min(pd[[i]]$fit, na.rm = TRUE)
            if (max(pd[[i]]$fit, na.rm = TRUE) > ylim[2]) 
              ylim[2] <- max(pd[[i]]$fit, na.rm = TRUE)
          }
        }
        if (partial.resids) {
          ul <- max(pd[[i]]$p.resid, na.rm = TRUE)
          if (ul > ylim[2]) 
            ylim[2] <- ul
          ll <- min(pd[[i]]$p.resid, na.rm = TRUE)
          if (ll < ylim[1]) 
            ylim[1] <- ll
        }
      }
  }
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
  if (n.para > 0) {
    class(x) <- c("gam", "glm", "lm")
    if (is.null(select)) {
      attr(x, "para.only") <- TRUE
      termplot(x, se = se, rug = rug, col.se = 1, col.term = 1, 
               main = attr(x$pterms, "term.labels"), ...)
    }
    else {
      if (select > m) {
        select <- select - m
        term.labels <- attr(x$pterms, "term.labels")
        term.labels <- term.labels[order == 1]
        if (select <= length(term.labels)) {
          termplot(x, terms = term.labels[select], se = se, 
                   rug = rug, col.se = 1, col.term = 1, ...)
        }
      }
    }
  }
  if (pages > 0) 
    par(oldpar)
}
<bytecode: 0x0000000013c6dc80>
  <environment: namespace:mgcv>