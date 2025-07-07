## ID: ccc.R, last updated 2025-07-01, F.Osorio

summary.L1ccc <- function(Nsize = 1000, type = "normal", eps = 0.05, vif = 0.10, alpha = 0.05) {
  stat <- matrix(0, nrow = 9, ncol = 6)
  SE   <- matrix(0, nrow = 9, ncol = 6)
  perc <- matrix(0, nrow = 9, ncol = 7)
  now <- proc.time()

  # m = 1
  rho <- 0.95
  cat(" 1/9:\n")
  set.seed(5)
  o <- simul.L1ccc(Nsize = Nsize, nobs = 25, type = type, rho = rho, eps = eps, vif = vif, alpha = alpha)
  perc[1,] <- o$percentage
  stat[1,] <- apply(o$stats, 2, mean)
  SE[1,] <- apply(o$SE, 2, mean)

  cat(" 2/9:\n")
  set.seed(23)
  o <- simul.L1ccc(Nsize = Nsize, nobs = 100, type = type, rho = rho, eps = eps, vif = vif, alpha = alpha)
  perc[2,] <- o$percentage
  stat[2,] <- apply(o$stats, 2, mean)
  SE[2,] <- apply(o$SE, 2, mean)

  cat(" 3/9:\n")
  set.seed(53)
  o <- simul.L1ccc(Nsize = Nsize, nobs = 400, type = type, rho = rho, eps = eps, vif = vif, alpha = alpha)
  perc[3,] <- o$percentage
  stat[3,] <- apply(o$stats, 2, mean)
  SE[3,] <- apply(o$SE, 2, mean)

  # m = 2
  rho <- 0.85
  cat(" 4/9:\n")
  set.seed(167)
  o <- simul.L1ccc(Nsize = Nsize, nobs = 25, type = type, rho = rho, eps = eps, vif = vif, alpha = alpha)
  perc[4,] <- o$percentage
  stat[4,] <- apply(o$stats, 2, mean)
  SE[4,] <- apply(o$SE, 2, mean)

  cat(" 5/9:\n")
  set.seed(239)
  o <- simul.L1ccc(Nsize = Nsize, nobs = 100, type = type, rho = rho, eps = eps, vif = vif, alpha = alpha)
  perc[5,] <- o$percentage
  stat[5,] <- apply(o$stats, 2, mean)
  SE[5,] <- apply(o$SE, 2, mean)

  cat(" 6/9:\n")
  set.seed(347)
  o <- simul.L1ccc(Nsize = Nsize, nobs = 400, type = type, rho = rho, eps = eps, vif = vif, alpha = alpha)
  perc[6,] <- o$percentage
  stat[6,] <- apply(o$stats, 2, mean)
  SE[6,] <- apply(o$SE, 2, mean)

  # m = 3
  rho <- 0.75
  cat(" 7/9:\n")
  set.seed(433)
  o <- simul.L1ccc(Nsize = Nsize, nobs = 25, type = type, rho = rho, eps = eps, vif = vif, alpha = alpha)
  perc[7,] <- o$percentage
  stat[7,] <- apply(o$stats, 2, mean)
  SE[7,] <- apply(o$SE, 2, mean)

  cat(" 8/9:\n")
  set.seed(577)
  o <- simul.L1ccc(Nsize = Nsize, nobs = 100, type = type, rho = rho, eps = eps, vif = vif, alpha = alpha)
  perc[8,] <- o$percentage
  stat[8,] <- apply(o$stats, 2, mean)
  SE[8,] <- apply(o$SE, 2, mean)

  cat(" 9/9:\n")
  set.seed(863)
  o <- simul.L1ccc(Nsize = Nsize, nobs = 400, type = type, rho = rho, eps = eps, vif = vif, alpha = alpha)
  perc[9,] <- o$percentage
  stat[9,] <- apply(o$stats, 2, mean)
  SE[9,] <- apply(o$SE, 2, mean)

  colnames(perc) <- c("Wald.N","Rao.N","Gradient.N","Wald.L","Rao.L","Gradient.L","T2")
  rownames(perc) <- c("1: 25","1:100","1:400","2: 25","2:100","2:400","3: 25","3:100","3:400")
  colnames(stat) <- c("ccc.N","L1.N","U.N","ccc.L","L1.L","U.L")
  rownames(stat) <- c("1: 25","1:100","1:400","2: 25","2:100","2:400","3: 25","3:100","3:400")
  colnames(SE) <- c("ccc.N","L1.N","U.N","ccc.L","L1.L","U.L")
  rownames(SE) <- c("1: 25","1:100","1:400","2: 25","2:100","2:400","3: 25","3:100","3:400")
  speed <- proc.time() - now

  list(stat = stat, SE = SE, perc = perc, speed  = speed) 
}

simul.L1ccc <- function(Nsize = 1000, nobs = 100, type = "normal", rho, eps = 0.05, vif = 0.10, alpha = 0.05)
{ ## function to perform the simulation experiment (Section 3 of the manuscript)
  stats <- matrix(0, nrow = Nsize, ncol = 6) # results container
  SE <- matrix(0, nrow = Nsize, ncol = 6) # results container
  tests <- matrix(0, nrow = Nsize, ncol = 7) # results container
  ok <- matrix(FALSE, nrow = Nsize, ncol = 7) # results container

  # setting parameters
  mu    <- c(0,0)
  Sigma <- matrix(c(1,rho,rho,1), ncol = 2)
  cutoff <- qchisq(1 - alpha, df = 1)

  pb <- txtProgressBar(min = 0, max = Nsize, style = 3)
  now <- proc.time()
  # Monte Carlo iterations
  for (i in 1:Nsize) {
    x <- switch(type,
                "normal" = rmnorm(n = nobs, mean = mu, Sigma = Sigma),
                "Laplace" = rmLaplace(n = nobs, center = mu, Scatter = Sigma),
                "Cauchy" = rmCauchy(n = nobs, center = mu, Scatter = Sigma),
                "CN" = rmcnorm(n = nobs, center = mu, Scatter = Sigma, epsilon = eps, vif = vif),
                stop("option not available."))
    # fitting postulated model and computing influence measures
    f0 <- ccc(x, method = "asymp", equal.means = TRUE)
    f1 <- l1ccc(x, equal.means = TRUE, boots = FALSE)

    # coefficient estimates
    stats[i,1] <- f0$Restricted$ccc
    stats[i,2] <- f0$Restricted$L1$rho1
    stats[i,3] <- f0$ustat$rhoc
    rho0 <- f1$Restricted$ccc
    stats[i,4] <- rho0
    stats[i,5] <- f1$Restricted$rho1
    stats[i,6] <- f1$ustat$rho1
    # standard errors
    SE[i,1] <- sqrt(f0$Restricted$var.ccc)
    SE[i,2] <- sqrt(f0$Restricted$L1$var.rho1)
    SE[i,3] <- sqrt(f0$ustat$var.rhoc)
    var.rho1 <- f1$Restricted$var.rho1
    SE[i,4] <- sqrt(var.rho1 * (1 - rho0))
    SE[i,5] <- sqrt(var.rho1)
    SE[i,6] <- sqrt(f1$ustat$var.rho1)
    # test statistics
    WN <- center.gaussian(f0, test = "Wald")$statistic
    RN <- center.gaussian(f0, test = "score")$statistic
    GN <- center.gaussian(f0, test = "gradient")$statistic
    WL <- center.test(f1, test = "Wald")$statistic
    RL <- center.test(f1, test = "score")$statistic
    GL <- center.test(f1, test = "gradient")$statistic
    T2 <- Hotelling(x)$statistic
    tests[i,1] <- WN
    tests[i,2] <- RN
    tests[i,3] <- GN
    tests[i,4] <- WL
    tests[i,5] <- RL
    tests[i,6] <- GL
    tests[i,7] <- T2
    # counts
    ok[i,1] <- WN > cutoff
    ok[i,2] <- RN > cutoff
    ok[i,3] <- GN > cutoff
    ok[i,4] <- WL > cutoff
    ok[i,5] <- RL > cutoff
    ok[i,6] <- GL > cutoff
    ok[i,7] <- T2 > cutoff

    # update progress bar
    setTxtProgressBar(pb, i)
  }
  close(pb)

  percentage <- rep(0, 7)
  percentage[1] <- sum(ok[,1], na.rm = TRUE) / Nsize
  percentage[2] <- sum(ok[,2], na.rm = TRUE) / Nsize
  percentage[3] <- sum(ok[,3], na.rm = TRUE) / Nsize
  percentage[4] <- sum(ok[,4], na.rm = TRUE) / Nsize
  percentage[5] <- sum(ok[,5], na.rm = TRUE) / Nsize
  percentage[6] <- sum(ok[,6], na.rm = TRUE) / Nsize
  percentage[7] <- sum(ok[,7], na.rm = TRUE) / Nsize

  colnames(stats) <- c("ccc.N","L1.N","U.N","ccc.L","L1.L","U.L")
  colnames(SE) <- c("ccc.N","L1.N","U.N","ccc.L","L1.L","U.L")
  colnames(tests) <- c("Wald.N","Rao.N","Gradient.N","Wald.L","Rao.L","Gradient.L","T2")
  colnames(ok) <- c("Wald.N","Rao.N","Gradient.N","Wald.L","Rao.L","Gradient.L","T2")
  names(percentage) <- c("Wald.N","Rao.N","Gradient.N","Wald.L","Rao.L","Gradient.L","T2")
  speed <- proc.time() - now
  out <- list(stats = stats, SE = SE, tests = tests, percentage = 100 * percentage, cutoff = cutoff, speed = speed)
  out
}
