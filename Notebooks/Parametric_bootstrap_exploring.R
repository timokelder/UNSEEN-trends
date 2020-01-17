##Boot doesnt work here
GEV_type = 'GEV'
extremes = extremes_wc
covariate = c(df_WC$Year) #1:35 instead of 1981:2015 for numerical efficiency
return_period = rperiods
covariate_values = c(1,35)
fit_Seas5_biascorr_nonstat <- fevd(extremes, type = GEV_type, location.fun = ~ covariate, ##Fitting the gev with a location and scale parameter linearly correlated to the covariate (years)
            scale.fun = ~ covariate, use.phi = TRUE)

params_matrix <- make.qcov(fit_Seas5_biascorr_nonstat, vals = list(mu1 = covariate_values,phi1 = covariate_values)) #Create a parameter matrix for the GEV fit
rvs=ci.fevd(fit,alpha = 0.05,type='return.level',return.period = return_period,method ="normal",qcov=params_matrix)  #Calculate the return values and confidence intervals for each year   

ci.fevd(fit_Seas5,alpha = 0.05,type='return.level',return.period = rperiods,method ="boot")

x=fit_Seas5
x=fit_obs
fit_Seas5_biascorr
const=is.fixedfevd(x)
qcov_1981 <- params_matrix[1,]
qcov_2015 <- params_matrix[2,]

R = 502 #Default number of bootstrap iterations 

x=fit_Seas5_biascorr_nonstat
Z <- rextRemes(x, n = R * x$n,qcov=params_matrix) #Draw samples from the extreme value distribution
Z_1981 <- matrix(Z[1,], x$n, R) # Create a matrix of a boorstrap per column
Z_2015 <- matrix(Z[2,], x$n, R) # Create a matrix of a boorstrap per column

Z <- rextRemes(x, n = R * x$n) #Draw samples from the extreme value distribution
Z <- matrix(Z, x$n, R) # Create a matrix of a boorstrap per column

###Fit the evd to each bootstrapped series
bfun <- function(z, x, y, p, ipar, eid, rate) {
  pm <- x$par.models
  if (is.null(y)) 
    fit <- fevd(x = z, threshold = x$threshold, location.fun = pm$location, 
                scale.fun = pm$scale, shape.fun = pm$shape, 
                use.phi = pm$log.scale, type = x$type, method = x$method, 
                initial = ipar, span = x$span, time.units = x$time.units, 
                period.basis = x$period.basis, optim.args = x$optim.args, 
                priorFun = x$priorFun, priorParams = x$priorParams, 
                verbose = FALSE)
  else fit <- fevd(x = z, data = y, threshold = x$threshold, 
                   location.fun = pm$location, scale.fun = pm$scale, 
                   shape.fun = pm$shape, use.phi = pm$log.scale, 
                   type = x$type, method = x$method, initial = ipar, 
                   span = x$span, time.units = x$time.units, period.basis = x$period.basis, 
                   optim.args = x$optim.args, priorFun = x$priorFun, 
                   priorParams = x$priorParams, verbose = FALSE)
  fit$cov.data <- y
  res <- distill(fit, cov = FALSE)
  return(res)
}

y <- datagrabber(x, response = FALSE)

pars <- apply(Z, 2, bfun, x = x, y = y, ipar = ipar)[1:np, 
                                                     


#rextremes: Understanding the sampling
p <- x$results$p
type <- x$type
n <- x$n
loc <- p["location"]
scale <- p["scale"]
shape <- p["shape"]
u <- 0

revd(n = n, loc = loc, scale = scale, shape = shape, 
     threshold = u, type = type)

#revd
z <- rexp(n) #Random generation for exponential distribution with rate =1
loc+scale*(z^(-shape) - 1)/shape ##what is z? return period? 
quantile(z,c(0.5,0.99))## it cant be return period

plot(1:100,1/1:100)

rep(loc, n) + sc[!id] * (z[!id]^(-sh[!id]) - 
                        1)/sh[!id]


                                                     ]
th.est <- numeric(3)
loc <- pars["location", ]
th.est[1] <- theta.hat["location"]

scale <- pars["scale", ]
th.est[2] <- theta.hat["scale"]
shape <- pars["shape", ]
th.est[3] <- theta.hat["shape"]
th <- rbind(loc, scale, shape)
rlfun <- function(theta, p, u, type, npy, rate) rlevd(period = p, 
                                                      loc = theta[1], scale = theta[2], shape = theta[3], 
                                                      threshold = u, type = type, npy = npy, rate = rate)
mod <- x$type
sam <- apply(th, 2, rlfun, p = return.period, u = x$threshold, 
             type = mod, npy = x$npy, rate = lam)


if (is.matrix(sam)) 
  rownames(sam) <- paste(rownames(sam), "-", 
                         x$period.basis, sep = "")
else sammy.name <- paste(return.period, "-", 
                         x$period.basis, sep = "")
theta.hat <- rlevd(period = return.period, loc = th.est[1], 
                   scale = th.est[2], shape = th.est[3], threshold = x$threshold, 
                   type = x$type, npy = x$npy, rate = lam)

if (is.matrix(sam)) {
  out <- apply(sam, 1, quantile, probs = c(alpha/2, 
                                           1 - alpha/2))
  out.names <- rownames(out)
  out <- rbind(out[1, ], theta.hat, out[2, ])
  rownames(out) <- c(out.names[1], "Estimate", 
                     out.names[2])
  colnames(out) <- rownames(sam)
  out <- t(out)
  attr(out, "data.name") <- x$call
  attr(out, "method") <- method.name
  attr(out, "conf.level") <- (1 - alpha) * 100
  attr(out, "R") <- R
  class(out) <- "ci"
  return(out)
}
else {
  out <- quantile(sam, probs = c(alpha/2, 1 - alpha/2))
  out <- c(out[1], mean(sam), out[2])
  attr(out, "R") <- R
}





###rlevd obtaining the return levels

theta.hat <- rlevd(period = return.period, loc = th.est[1], 
                   scale = th.est[2], shape = th.est[3], threshold = x$threshold, 
                   type = x$type, npy = x$npy, rate = lam)


function (n, loc = 0, scale = 1, shape = 0, threshold = 0, type = c("GEV", 
                                                                    "GP")) 
{
  type <- match.arg(type)
  type <- tolower(type)
  if (type == "gev") 
    z <- rexp(n)
  else if (type == "gp") {
    z <- runif(n)
    loc <- threshold
  }
  else stop("revd: invalid type argument.")
  out <- numeric(n) + NA
  if (min(scale) < 0) 
    stop("revd: scale parameter(s) must be positively valued.")
  if (length(loc) == 1) 
    loc <- rep(loc, n)
  if (length(scale) == 1) 
    sc <- rep(scale, n)
  else sc <- scale
  if (length(shape) == 1) 
    sh <- rep(shape, n)
  else sh <- shape
  if (length(loc) != n || length(sc) != n || length(sh) != 
      n) 
    stop("revd: parameters must have length equal to 1 or n.")
  id <- sh == 0
  if (any(id)) {
    if (type == "gev") 
      out[id] <- loc[id] - sc[id] * log(z[id])
    else if (type == "gp") 
      out[id] <- loc[id] + rexp(n, rate = 1/sc[id])
  }
  if (any(!id)) 
    out[!id] <- loc[!id] + sc[!id] * (z[!id]^(-sh[!id]) - 
                                        1)/sh[!id]
  return(out)
}


  type <- x$type
  if (is.fixedfevd(x)) {
    if (missing(n)) 
      n <- x$n
    p <- x$results$p
    pnames <- names(p)
    if (is.element("log.scale", pnames)) {
      id <- pnames == "log.scale"
      p[id] <- exp(p[id])
      pnames[id] <- "scale"
      names(p) <- pnames
    }
    scale <- p["scale"]
    if (x$par.models$log.scale) 
      scale <- exp(scale)
    if (type == "GEV") {
      loc <- p["location"]
      shape <- p["shape"]
      u <- 0
    }
    else if (type == "GP") {
      loc <- 0
      u <- x$threshold
      shape <- p["shape"]
    }
    else if (type == "PP") {
      loc <- p["location"]
      u <- x$threshold
      shape <- p["shape"]
      scale <- scale + shape * (u - loc)
      type <- "GP"
    }
    else if (type == "Gumbel") {
      loc <- p["location"]
      shape <- 0
      type <- "GEV"
      u <- 0
    }
    else if (is.element(type, c("Weibull", "Frechet"))) {
      loc <- p["location"]
      shape <- p["shape"]
      u <- 0
      type <- "GEV"
    }
    else if (type == "Exponential") {
      loc <- 0
      shape <- 0
      u <- x$threshold
      type <- "GP"
    }
    else if (is.element(type, c("Beta", "Pareto"))) {
      loc <- 0
      shape <- p["shape"]
      u <- x$threshold
      type <- "GP"
    }
    return(revd(n = n, loc = loc, scale = scale, shape = shape, 
                threshold = u, type = type))
  }
  else {
    if (missing(n)) 
      n <- 1
    p <- findpars(x, qcov = qcov)
    if (is.null(qcov)) 
      K <- x$n
    else K <- dim(qcov)[1]
    if (!is.null(qcov)) 
      u <- qcov[, "threshold"]
    else if (is.null(x$threshold)) 
      u <- rep(0, K)
    else u <- x$threshold
    if (is.null(p$location)) 
      loc <- rep(0, K)
    else loc <- p$location
    if (is.null(p$shape)) {
      if (type == "Gumbel") 
        shape <- rep(0, K)
      else shape <- rep(1e-10, K)
    }
    else shape <- p$shape
    scale <- p$scale
    if (type == "PP") 
      scale <- scale + shape * (u - loc)
    theta <- cbind(u, loc, scale, shape)
    if (is.null(qcov)) 
      y <- c(datagrabber(x, cov.data = FALSE))
    if (is.element(type, c("GEV", "Gumbel", "Weibull", 
                           "Frechet"))) {
      out <- matrix(NA, K, n)
      for (i in 1:n) {
        z <- revd(K, type = "GEV")
        out[, i] <- revtrans.evd(z = z, threshold = theta[, 
                                                          1], location = theta[, 2], scale = theta[, 
                                                                                                   3], shape = theta[, 4], type = type)
      }
    }
    else {
      if (is.null(qcov)) {
        eid <- y > u
        m <- sum(eid)
        theta <- theta[eid, ]
      }
      else {
        m <- K
      }
      out <- matrix(NA, m, n)
      if (type == "PP") 
        type2 <- "GP"
      else type2 <- type
      for (i in 1:n) {
        z <- revd(m, type = "GP")
        out[, i] <- revtrans.evd(z = z, threshold = theta[, 
                                                          1], location = theta[, 2], scale = theta[, 
                                                                                                   3], shape = theta[, 4], type = type2)
      }
    }
    return(out)
  }
})

ci.fevd.mle
function (x, alpha = 0.05, type = c("return.level", "parameter"), 
          return.period = 100, which.par = 1, R = 502, method = c("normal", 
                                                                  "boot", "proflik"), xrange = NULL, nint = 20, 
          verbose = FALSE, tscale = FALSE, return.samples = FALSE, 
          ...) 
{
  if (missing(method)) 
    miss.meth <- TRUE
  else miss.meth <- FALSE
  method <- tolower(method)
  method <- match.arg(method)
  type <- tolower(type)
  type <- match.arg(type)
  theta.hat <- x$results$par
  theta.names <- names(theta.hat)
  np <- length(theta.hat)
  if (type == "parameter" && missing(which.par)) 
    which.par <- 1:np
  if (any(theta.names == "log.scale")) {
    id <- theta.names == "log.scale"
    theta.hat[id] <- exp(theta.hat[id])
    theta.names[id] <- "scale"
    names(theta.hat) <- theta.names
  }
  const <- is.fixedfevd(x)
  if (type == "return.level") 
    par.name <- paste(return.period, "-", x$period.basis, 
                      " return level", sep = "")
  else if (type == "parameter") 
    par.name <- theta.names[which.par]
  if (type == "return.level" && !const) {
    return(ci.rl.ns.fevd.mle(x = x, alpha = alpha, return.period = return.period, 
                             method = method, verbose = verbose, return.samples = return.samples, 
                             ...))
  }
  if (type == "parameter") 
    p <- theta.hat[which.par]
  else {
    if (is.element(x$type, c("PP", "GP", "Beta", 
                             "Pareto", "Exponential"))) 
      lam <- mean(c(datagrabber(x)[, 1]) > x$threshold)
    else lam <- 1
    if (is.element(x$type, c("PP", "GEV", "Gumbel", 
                             "Weibull", "Frechet"))) 
      loc <- theta.hat["location"]
    else loc <- 0
    scale <- theta.hat["scale"]
    if (!is.element(x$type, c("Gumbel", "Exponential"))) 
      shape <- theta.hat["shape"]
    else shape <- 0
    if (x$type == "PP") 
      mod <- "GEV"
    else mod <- x$type
    p <- rlevd(period = return.period, loc = loc, scale = scale, 
               shape = shape, threshold = x$threshold, type = mod, 
               npy = x$npy, rate = lam)
  }
  if (verbose) {
    cat("\n", "Preparing to calculate ", (1 - 
                                            alpha) * 100, "% CI for ", ifelse(type == "return.level", 
                                                                              paste(return.period, "-", x$period.basis, " return level", 
                                                                                    sep = ""), paste(par.name, " parameter", 
                                                                                                     sep = "")), "\n")
    cat("\n", "Model is ", ifelse(const, " fixed", 
                                  " non-stationary."), "\n")
    if (method == "normal") 
      cat("\n", "Using Normal Approximation Method.\n")
    else if (method == "boot") 
      cat("\n", "Using Bootstrap Method.\n")
    else if (method == "proflik") 
      cat("\n", "Using Profile Likelihood Method.\n")
  }
  if (method == "normal") {
    method.name <- "Normal Approx."
    z.alpha <- qnorm(alpha/2, lower.tail = FALSE)
    cov.theta <- parcov.fevd(x)
    if (is.null(cov.theta)) 
      stop("ci: Sorry, unable to calculate the parameter covariance matrix.  Maybe try a different method.")
    var.theta <- diag(cov.theta)
    if (any(var.theta < 0)) 
      stop("ci: negative Std. Err. estimates obtained.  Not trusting any of them.")
    if (type == "parameter") {
      se.theta <- sqrt(var.theta)
      if (tscale) {
        if (!const && !is.element("scale", theta.names) && 
            !is.element("shape", theta.names) && 
            !all(x$threshold == x$threshold[1])) {
          stop("ci: invalid argument configurations.")
        }
        if (!is.element(x$type, c("GP", "Beta", 
                                  "Pareto"))) 
          stop("ci: invalid argument configurations.")
        theta.hat["scale"] <- theta.hat["scale"] - 
          theta.hat["shape"] * x$threshold
        theta.names[theta.names == "scale"] <- "tscale"
        if (!any(theta.names[which.par] == "tscale")) 
          stop("ci: invalid argument configurations.")
        names(theta.hat) <- theta.names
        p <- theta.hat[which.par]
        d <- rbind(1, -x$threshold)
        names(se.theta) <- theta.names
        se.theta["tscale"] <- sqrt(t(d) %*% cov.theta %*% 
                                     d)
      }
      else se.theta <- sqrt(var.theta)[which.par]
      se.theta <- se.theta[which.par]
      par.name <- theta.names[which.par]
    }
    else if (type == "return.level") {
      grads <- rlgrad.fevd(x, period = return.period)
      grads <- t(grads)
      if (is.element(x$type, c("GP", "Beta", 
                               "Pareto", "Exponential"))) {
        if (x$type == "Exponential") 
          cov.theta <- diag(c(lam * (1 - lam)/x$n, var.theta))
        else cov.theta <- rbind(c(lam * (1 - lam)/x$n, 
                                  0, 0), cbind(0, cov.theta))
      }
      else lam <- 1
      var.theta <- t(grads) %*% cov.theta %*% grads
    }
    else stop("ci: invalid type argument.  Must be return.level or parameter.")
    if (length(p) > 1) {
      if (type == "return.level") 
        se.theta <- sqrt(diag(var.theta))
      out <- cbind(p - z.alpha * se.theta, p, p + z.alpha * 
                     se.theta)
      rownames(out) <- par.name
      conf.level <- paste(round((1 - alpha) * 100, digits = 2), 
                          "%", sep = "")
      colnames(out) <- c(paste(conf.level, " lower CI", 
                               sep = ""), "Estimate", paste(conf.level, 
                                                            " upper CI", sep = ""))
      attr(out, "data.name") <- x$call
      attr(out, "method") <- method.name
      attr(out, "conf.level") <- (1 - alpha) * 100
      class(out) <- "ci"
      return(out)
    }
    else out <- c(p - z.alpha * sqrt(var.theta[which.par]), 
                  p, p + z.alpha * sqrt(var.theta[which.par]))
  }
  else if (method == "boot") {
    method.name <- "Parametric Bootstrap"
    if (verbose) 
      cat("\n", "Simulating data from fitted model.  Size = ", 
          R, "\n")
    if (const) {
      if (is.null(x$blocks)) {
        Z <- rextRemes(x, n = R * x$n)
        Z <- matrix(Z, x$n, R)
      }
      else {
        Z <- rextRemes(x, n = round(R * x$npy * x$blocks$nBlocks))
        Z <- matrix(Z, round(x$npy * x$blocks$nBlocks), 
                    R)
      }
    }
    else Z <- rextRemes(x, n = R)
    if (verbose) 
      cat("\n", "Simulated data found.\n")
    y <- datagrabber(x, response = FALSE)
    if (is.element(x$type, c("PP", "GP", "Exponential", 
                             "Beta", "Pareto"))) {
      x2 <- datagrabber(x, cov.data = FALSE)
      eid <- x2 > x$threshold
      Z2 <- matrix(x$threshold, x$n, R)
      Z2[eid, ] <- Z[eid, ]
      Z <- Z2
      lam <- mean(eid)
    }
    else {
      eid <- !logical(x$n)
      lam <- 1
    }
    ipar <- list()
    if (any(is.element(c("location", "mu0"), 
                       theta.names))) {
      if (is.element("location", theta.names)) 
        ipar$location <- theta.hat["location"]
      else {
        id <- substring(theta.names, 1, 2) == "mu"
        ipar$location <- theta.hat[id]
      }
    }
    if (is.element("scale", theta.names)) 
      ipar$scale <- theta.hat["scale"]
    else {
      if (!x$par.models$log.scale) 
        id <- substring(theta.names, 1, 3) == "sig"
      else id <- substring(theta.names, 1, 3) == "phi"
      ipar$scale <- theta.hat[id]
    }
    if (!is.element(x$type, c("Gumbel", "Exponential"))) {
      if (is.element("shape", theta.names)) 
        ipar$shape <- theta.hat["shape"]
      else {
        id <- substring(theta.names, 1, 2) == "xi"
        ipar$shape <- theta.hat[id]
      }
    }
    bfun <- function(z, x, y, p, ipar, eid, rate) {
      pm <- x$par.models
      if (is.null(y)) 
        fit <- fevd(x = z, threshold = x$threshold, location.fun = pm$location, 
                    scale.fun = pm$scale, shape.fun = pm$shape, 
                    use.phi = pm$log.scale, type = x$type, method = x$method, 
                    initial = ipar, span = x$span, time.units = x$time.units, 
                    period.basis = x$period.basis, optim.args = x$optim.args, 
                    priorFun = x$priorFun, priorParams = x$priorParams, 
                    verbose = FALSE)
      else fit <- fevd(x = z, data = y, threshold = x$threshold, 
                       location.fun = pm$location, scale.fun = pm$scale, 
                       shape.fun = pm$shape, use.phi = pm$log.scale, 
                       type = x$type, method = x$method, initial = ipar, 
                       span = x$span, time.units = x$time.units, period.basis = x$period.basis, 
                       optim.args = x$optim.args, priorFun = x$priorFun, 
                       priorParams = x$priorParams, verbose = FALSE)
      fit$cov.data <- y
      res <- distill(fit, cov = FALSE)
      return(res)
    }
    if (verbose) 
      cat("\n", "Fitting model to simulated data sets (this may take a while!).")
    if (type == "parameter") {
      sam <- apply(Z, 2, bfun, x = x, y = y, ipar = ipar)
      if (tscale) {
        if (!const && !is.element("scale", theta.names) && 
            !is.element("shape", theta.names)) 
          stop("ci: invalid argument configurations.")
        if (!is.element(x$type, c("GP", "Beta", 
                                  "Pareto"))) 
          stop("ci: invalid argument configurations.")
        sam["scale", ] <- sam["scale", ] - 
          sam["shape", ] * x$threshold
        theta.hat["scale"] <- theta.hat["scale"] - 
          theta.hat["shape"] * x$threshold
        theta.names[theta.names == "scale"] <- "tscale"
        rownames(sam) <- theta.names
        names(theta.hat) <- theta.names
      }
      sam <- sam[which.par, ]
      if (return.samples) 
        return(t(sam))
    }
    else if (type == "return.level") {
      pars <- apply(Z, 2, bfun, x = x, y = y, ipar = ipar)[1:np, 
                                                           ]
      th.est <- numeric(3)
      if (is.element(x$type, c("PP", "GEV", 
                               "Gumbel", "Weibull", "Frechet"))) {
        loc <- pars["location", ]
        th.est[1] <- theta.hat["location"]
      }
      else loc <- rep(0, R)
      scale <- pars["scale", ]
      th.est[2] <- theta.hat["scale"]
      if (!is.element(x$type, c("Gumbel", "Exponential"))) {
        shape <- pars["shape", ]
        th.est[3] <- theta.hat["shape"]
      }
      else {
        shape <- rep(1e-10, R)
        th.est[3] <- 1e-08
      }
      if (return.samples) 
        out <- t(pars)
      th <- rbind(loc, scale, shape)
      rlfun <- function(theta, p, u, type, npy, rate) rlevd(period = p, 
                                                            loc = theta[1], scale = theta[2], shape = theta[3], 
                                                            threshold = u, type = type, npy = npy, rate = rate)
      if (x$type == "PP") 
        mod <- "GEV"
      else mod <- x$type
      sam <- apply(th, 2, rlfun, p = return.period, u = x$threshold, 
                   type = mod, npy = x$npy, rate = lam)
      if (is.matrix(sam)) 
        rownames(sam) <- paste(rownames(sam), "-", 
                               x$period.basis, sep = "")
      else sammy.name <- paste(return.period, "-", 
                               x$period.basis, sep = "")
      if (return.samples) {
        if (is.matrix(sam)) 
          out <- cbind(pars, t(sam))
        else {
          onames <- colnames(out)
          out <- cbind(out, sam)
          colnames(out) <- c(onames, sammy.name)
        }
        return(out)
      }
      theta.hat <- rlevd(period = return.period, loc = th.est[1], 
                         scale = th.est[2], shape = th.est[3], threshold = x$threshold, 
                         type = x$type, npy = x$npy, rate = lam)
    }
    else stop("ci: invalid type argument.  Must be return.level or parameter.")
    if (is.matrix(sam)) {
      out <- apply(sam, 1, quantile, probs = c(alpha/2, 
                                               1 - alpha/2))
      out.names <- rownames(out)
      out <- rbind(out[1, ], theta.hat, out[2, ])
      rownames(out) <- c(out.names[1], "Estimate", 
                         out.names[2])
      colnames(out) <- rownames(sam)
      out <- t(out)
      attr(out, "data.name") <- x$call
      attr(out, "method") <- method.name
      attr(out, "conf.level") <- (1 - alpha) * 100
      attr(out, "R") <- R
      class(out) <- "ci"
      return(out)
    }
    else {
      out <- quantile(sam, probs = c(alpha/2, 1 - alpha/2))
      out <- c(out[1], mean(sam), out[2])
      attr(out, "R") <- R
    }
    if (verbose) 
      cat("\n", "Finished fitting model to simulated data.\n")
  }
  else if (method == "proflik") {
    if (x$type == "PP" && !is.null(x$blocks)) 
      stop("ci: cannot do profile likelihood with blocks.")
    if (tscale) 
      stop("ci: invalid argument configurations.")
    if (type == "parameter" && length(which.par) > 
        1) 
      stop("ci: can only do one parameter at a time with profile likelihood method.")
    else if (type == "return.level" && length(return.period) > 
             1) 
      stop("ci: can only do one return level at a time with profile likelihood method.")
    method.name <- "Profile Likelihood"
    if (verbose) {
      if (x$type != "PP") 
        cat("\n", "Calculating profile likelihood.  This may take a few moments.\n")
      else cat("\n", "Calculating profile likelihood.  This may take several moments.\n")
    }
    if (is.null(xrange)) {
      hold2 <- c(ci(x, alpha = alpha, method = "normal", 
                    type = type, return.period = return.period, which.par = which.par))[c(1, 
                                                                                          3)]
      if (!any(is.na(hold2))) 
        xrange <- range(c(hold2, log2(hold2), 4 * hold2, 
                          hold2 - 4 * hold2, hold2 + 4 * hold2), finite = TRUE)
      else if (!is.na(hold2[2])) 
        xrange <- range(c(p - 2 * abs(log2(abs(p))), 
                          hold2[2], 4 * hold2[2], -4 * hold2[2], log2(p)), 
                        finite = TRUE)
      else if (!is.na(hold2[1])) 
        xrange <- range(c(p - 2 * abs(log2(abs(p))), 
                          hold2[1], 4 * hold2[1], -4 * hold2[1], log2(p)), 
                        finite = TRUE)
      else if (all(is.na(hold2))) 
        xrange <- c(p - 2 * abs(log2(abs(p))), p + 2 * 
                      abs(log2(abs(p))))
      if (verbose) 
        cat("\n", "Using a range of ", xrange[1], 
            " to ", xrange[2], "\n")
    }
    if (is.null(x$blocks)) {
      if (!is.null(xrange)) 
        hold <- profliker(x, type = type, xrange = xrange, 
                          return.period = return.period, which.par = which.par, 
                          nint = nint, plot = verbose, ...)
      else hold <- profliker(x, type = type, return.period = return.period, 
                             which.par = which.par, nint = nint, plot = verbose, 
                             ...)
    }
    else stop("Sorry: profile likelihood with blocks is not supported.")
    ma <- -x$results$value
    crit <- ma - 0.5 * qchisq(1 - alpha, 1)
    if (verbose) {
      cat("\n", "Profile likelihood has been calculated.  Now, trying to find where it crosses the critical value = ", 
          crit, "\n")
      abline(h = crit, col = "blue")
    }
    crit2 <- ma - 0.5 * qchisq((1 - alpha) + abs(log2(1 - 
                                                        alpha))/2, 1)
    id <- hold > crit2
    z <- seq(xrange[1], xrange[2], length = length(hold))
    z <- z[id]
    parlik <- hold[id]
    smth <- spline(z, parlik, n = 200)
    ind <- smth$y > crit
    out <- range(smth$x[ind])
    if (verbose) 
      abline(v = out, lty = 2, col = "darkblue", 
             lwd = 2)
    out <- c(out[1], p, out[2])
  }
  else stop("ci: invalid method argument.")
  conf.level <- paste(round((1 - alpha) * 100, digits = 2), 
                      "%", sep = "")
  names(out) <- c(paste(conf.level, " lower CI", sep = ""), 
                  par.name, paste(conf.level, " upper CI", sep = ""))
  attr(out, "data.name") <- x$call
  attr(out, "method") <- method.name
  attr(out, "conf.level") <- (1 - alpha) * 100
  class(out) <- "ci"
  return(out)
}