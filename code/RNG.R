## ID: RNG.R, last updated 2025-06-26, F.Osorio

rmCauchy <- 
function(n = 1, center = rep(0, nrow(Scatter)), Scatter = diag(length(center)))
{ # multivariate Cauchy random generation
  if (n <= 0)
    stop("n must be a positive integer")
  if (length(center) != nrow(Scatter))
    stop("center and scatter have non-conforming size")
  p <- nrow(Scatter)

  y <- matrix(0, nrow = n, ncol = p)
  # call C code
  y <- .C("RNG_cauchy", 
          y = as.double(y),
          n = as.integer(n),
          p = as.integer(p),
          center = as.double(center),
          Scatter = as.double(Scatter))$y
  y <- matrix(y, nrow = n, byrow = TRUE)
  y
}

rmcnorm <- 
function(n = 1, center = rep(0, nrow(Scatter)), Scatter = diag(length(center)), 
  epsilon = 0.05, vif = 0.25)
{
  if (n <= 0)
    stop("n must be a positive integer")
  if (length(center) != nrow(Scatter))
    stop("center and scatter have non-conforming size")
  if ((epsilon < 0) || (epsilon > 1))
      stop("contamination percentage must be in [0,1]")
  if ((vif <= 0) || (vif >= 1))
      stop("variance inflation factor must be in (0,1)")
  p <- nrow(Scatter)

  y <- matrix(0, nrow = n, ncol = p)
  dy <- dim(y)
  # call C code
  y <- .C("RNG_contaminated", 
          y = as.double(y),
          n = as.integer(dy[1]),
          p = as.integer(dy[2]),
          center = as.double(center),
          Scatter = as.double(Scatter),
          epsilon = as.double(epsilon),
          vif = as.double(vif))$y
  y <- matrix(y, nrow = n, byrow = TRUE)
  y
}
