## ID: envelope.R, last updated 2026-01-02, F.Osorio

qenvelope <- function(object, reps = 50, conf = 0.95, plot.it = TRUE)
{ ## simulated envelope based on 'quantile residuals'
  ok <- inherits(object, "LaplaceFit") || inherits(object, "studentFit")
  if (!ok)
    stop("Use only with 'studentFit' or 'LaplaceFit' objects")
  if (inherits(object, "studentFit")) {
    zero.eta <- (object$eta == 0)
    if (!zero.eta)
      stop("Only for use with 'eta = 0'")
  }

  qtrans.gaussian <- function(distances, p) {
    distances <- sqrt(distances)
    z <- qnorm(pchi(distances, df = p))
    z
  }

  qtrans.Laplace <- function(distances, p) {
    z <- qnorm(pgamma(distances, shape = p, scale = 2.0))
    z
  }
  
  envel <- function(n, mean, Sigma, type, reps, conf) {
    conf <- 1 - conf
    # initialize progress bar
    cat(" Progress:\n")
    pb <- txtProgressBar(min = 0, max = reps, style = 3)
    elims <- matrix(0, nrow = n, ncol = reps)
    for (i in 1:reps) {
      switch(type, 
             "normal" = {
              x <- rmt(n, mean = mean, Sigma = Sigma, eta = 0)
              fit <- studentFit(x, family = Student(eta = 0))
              z <- qtrans.gaussian(fit$distances, p)
             },
             "Laplace" = {
              x <- rmLaplace(n, center = mean, Scatter = Sigma)
              fit <- LaplaceFit(x)
              z <- qtrans.Laplace(fit$distances, p)
             })
      elims[,i] <- sort(z)
      # update progress bar
      setTxtProgressBar(pb, i)
    }
    close(pb)
    band <- matrix(0, nrow = n, ncol = 2)
    for (i in 1:n)
      band[i,] <- quantile(elims[i,], probs = c(conf / 2, 1 - conf / 2))
    band
  }

  if (inherits(object, "studentFit"))
    dist <- "normal"
  if (inherits(object, "LaplaceFit"))
    dist <- "Laplace"
 
  n <- object$dims[1]
  p <- object$dims[2]
  z <- switch(dist, 
              "normal" = qtrans.gaussian(object$distances, p),
              "Laplace" = qtrans.Laplace(object$distances, p))

  if (plot.it) {
    band  <- envel(n, object$center, object$Scatter, dist, reps, conf)
    ylim <- range(z, band)
    qqnorm(z, ylim = ylim, main = "Quantile distances Q-Q plot")
    par(new = TRUE)
    qqnorm(band[,1], axes = F, main = "", xlab = "", ylab = "", ylim = ylim, type = "l", lwd = 2, col = "red")
    par(new = TRUE)
    qqnorm(band[,2], axes = F, main = "", xlab = "", ylab = "", ylim = ylim, type = "l", lwd = 2, col = "red")
  }

  output <- list(transformed = z)
  if (plot.it)
    output$envelope = band
  invisible(output)
}
