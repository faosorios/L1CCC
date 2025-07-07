## ID: center_test.R, last updated 2025-06-27, F.Osorio

center.test <-
function(o, test = "Wald")
{
  ## local functions
  restrictedScore <-
  function(z, dims, weights, center)
  { # compute the unscaled score
    n <- dims[1]
    p <- dims[2]
    y <- z - matrix(rep(center, times = n), ncol = p, byrow = TRUE)
    y <- weights * y
    s <- apply(y, 2, sum)
    s
  }

  # checking 
  if (!inherits(o, "L1ccc"))
    stop("Use only with 'L1ccc' objects")
  dx <- o$dims
  n <- dx[1]
  p <- dx[2]
  if (p != 2) stop("only implemented for p = 2.")

  # extracting elements
  x <- o$x
  center <- o$center
  Scatter <- o$Scatter
  mu0 <- o$Restricted$Fitted$center
  Sigma0 <- o$Restricted$Fitted$Scatter
  wts0 <- o$Restricted$Fitted$weights
  null.fit <- list(center = mu0, Scatter = Sigma0, weights = wts0)

  switch(test,
    Wald = {
      A <- matrix(c(1,-1), ncol = 2)
      delta <- c(A %*% center)
      s2 <- c(A %*% Scatter %*% t(A))
      stat <- (0.25 * n / p) * (delta^2) / s2
      names(stat) <- "Wald"
      method <- "Wald test"
    },
    score = {
      s <- restrictedScore(x, dx, wts0, mu0)
      stat <- (4 * p / n) * c(t(s) %*% solve(Sigma0, s))
      names(stat)<-"Score"
      method <- "Score test"
    },
    gradient = {
      delta <- solve(Sigma0, center)
      s <- restrictedScore(x, dx, wts0, mu0)
      stat <- sum(s * delta)
      names(stat) <- "Gradient"
      method <- "Gradient test"
    },
    stop(paste("unimplemented test:", test)))
  pval <- 1 - pchisq(stat, df = 1)
  
  ## output object
  delta <- center[1] - center[2]
  z <- list(statistic = stat, parameter = 1, p.value = pval, estimate = delta,
            null.value = 0, method = method, null.fit = null.fit)
  class(z) <- "center.test"
  z
}

center.gaussian <-
function(o, test = "Wald")
{
  ## local functions
  restrictedScore <-
  function(z, dims, center)
  { # compute the unscaled score
    n <- dims[1]
    p <- dims[2]
    y <- z - matrix(rep(center, times = n), ncol = p, byrow = TRUE)
    s <- apply(y, 2, sum)
    s
  }

  # checking 
  if (!inherits(o, "ccc"))
    stop("Use only with 'ccc' objects")
  dx <- o$dims
  n <- dx[1]
  p <- dx[2]
  if (p != 2) stop("only implemented for p = 2.")

  # extracting elements
  x <- o$x
  xbar <- o$center
  S <- o$cov
  mu0 <- o$Restricted$center
  Sigma0 <- o$Restricted$cov
  null.fit <- list(mean = mu0, cov = Sigma0)

  switch(test,
    Wald = {
      A <- matrix(c(1,-1), ncol = 2)
      delta <- c(A %*% xbar)
      s2 <- c(A %*% S %*% t(A))
      stat <- n * (delta^2) / s2
      names(stat) <- "Wald"
      method <- "Wald test"
    },
    score = {
      s <- restrictedScore(x, dx, mu0)
      stat <- (1 / n) * c(t(s) %*% solve(Sigma0, s))
      names(stat)<-"Score"
      method <- "Score test"
    },
    gradient = {
      delta <- solve(Sigma0, xbar)
      s <- restrictedScore(x, dx, mu0)
      stat <- sum(s * delta)
      names(stat) <- "Gradient"
      method <- "Gradient test"
    },
    stop(paste("unimplemented test:", test)))
  pval <- 1 - pchisq(stat, df = 1)
  
  ## output object
  delta <- xbar[1] - xbar[2]
  z <- list(statistic = stat, parameter = 1, p.value = pval, estimate = delta,
            null.value = 0, method = method, null.fit = null.fit)
  class(z) <- "center.test"
  z
}

print.center.test <- function(x, digits = 4, ...)
{
  cat("\n")
  cat(paste(x$method, "for the equality of location parameters", sep = " "), "\n")
  cat("\n")
  out <- character()
  out <- c(out, paste(names(x$statistic), "statistic =",
                      format(round(x$statistic, digits = digits))))
  out <- c(out, paste("df =", x$parameter))
  out <- c(out, paste("p-value =", format(round(x$p.value, digits = digits))))
  cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
  cat(paste("alternative hypothesis: true difference parameter is not equal to:", format(round(x$null.value, digits = digits)), sep = " "), "\n")
  cat(paste("sample estimate:", format(round(x$estimate, digits = digits)), sep = " "), "\n")
  invisible(x)
}

Hotelling <- function(x) 
{ ## generalized Hotelling's T-squared

  # checking 
  p <- ncol(x)
  if (p != 2) stop("only implemented for p = 2.")

  # computing sample mean and covariance
  xbar <- colMeans(x)
  S <- cov(x)

  # computing Hotelling statistic
  n <- nrow(x)
  A <- matrix(c(1,-1), ncol = 2)
  delta <- c(A %*% xbar)
  s2 <- c(A %*% S %*% t(A))
  stat <- n * (delta^2) / s2
  names(stat) <- "Hotelling"
  method <- "Hotelling's T-squared test"
  pval <- 1 - pchisq(stat, df = 1)
  
  ## output object
  z <- list(statistic = stat, parameter = 1, p.value = pval, estimate = delta,
            null.value = 0, method = method)
  class(z) <- "Hotelling"
  z
}

print.Hotelling <- function(x, digits = 4, ...)
{
  cat("\n")
  cat(paste(x$method, "for the equality of location parameters", sep = " "), "\n")
  cat("\n")
  out <- character()
  out <- c(out, paste(names(x$statistic), "statistic =",
                      format(round(x$statistic, digits = digits))))
  out <- c(out, paste("df =", x$parameter))
  out <- c(out, paste("p-value =", format(round(x$p.value, digits = digits))))
  cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
  cat(paste("alternative hypothesis: true difference parameter is not equal to:", format(round(x$null.value, digits = digits)), sep = " "), "\n")
  cat(paste("sample estimate:", format(round(x$estimate, digits = digits)), sep = " "), "\n")
  invisible(x)
}
