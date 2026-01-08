## script related to Supplementary material

## loading packages and reading sources
library(L1pack)
library(MVT)
source("../code/envelope.R")

## Fig. 1(a)
env0 <- qenvelope(fm0, reps = 5000)
ylim <- c(-3.3, 4.2)
par(pty = "s", mai = c(1,1,.35,.35))
qqnorm(env0$trans, ylim = ylim, main = "", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", lwd = 2, cex.lab = 1.3)
par(new = TRUE)
qqnorm(env0$envelope[,1], axes = FALSE, ylim = ylim, main = "", xlab = "", ylab = "", lwd = 2, col = "red", type = "l")
par(new = TRUE)
qqnorm(env0$envelope[,2], axes = FALSE, ylim = ylim, main = "", xlab = "", ylab = "", lwd = 2, col = "red", type = "l")

## Fig. 1(b)
env1 <- qenvelope(fm1, reps = 5000)
ylim <- c(-3.3, 4)
par(pty = "s", mai = c(1,1,.35,.35))
qqnorm(env1$trans, ylim = ylim, main = "", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", lwd = 2, cex.lab = 1.3)
par(new = TRUE)
qqnorm(env1$envelope[,1], axes = FALSE, ylim = ylim, main = "", xlab = "", ylab = "", lwd = 2, col = "red", type = "l")
par(new = TRUE)
qqnorm(env1$envelope[,2], axes = FALSE, ylim = ylim, main = "", xlab = "", ylab = "", lwd = 2, col = "red", type = "l")
