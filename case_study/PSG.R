## loading packages, reading sources and attaching dataset
library(fastmatrix)
library(L1pack)
library(MVT)
source("../code/center_test.R")
data(PSG)

## model fit and simulated envelope under normal distribution
fm0 <- studentFit(~ manual + automated, data = PSG, family = Student(eta = 0))
env <- envelope.student(fm0, reps = 5000)

## Fig. 1(a)
ylim <- c(-2.5, 4.15)
par(pty = "s", mai = c(1,1,.35,.35))
qqnorm(env$trans, ylim = ylim, main = "", lwd = 2, cex.lab = 1.3)
par(new = TRUE)
qqnorm(env$envelope[,1], axes = FALSE, ylim = ylim, main = "", xlab = "", ylab = "", lwd = 2, col = "red", type = "l")
par(new = TRUE)
qqnorm(env$envelope[,2], axes = FALSE, ylim = ylim, main = "", xlab = "", ylab = "", lwd = 2, col = "red", type = "l")

## Fig. 1(b)
z0 <- ccc(~ manual + automated, data = PSG, method = "asymp", equal.means = TRUE)
nobs <- nrow(PSG)

# removing i-th observation
rhoc <- rep(0, nobs)
for (i in 1:nobs)
  rhoc[i] <- ccc(~ manual + automated, data = PSG, subset = -i)$ccc

obs <- c(1,30,79)
cutoff <- z0$ccc
par(pty = "s", mai = c(1,1,.35,.35))
plot(rhoc, type = "b", ylim = c(.65, .755), ylab = "CCC estimate", lwd = 2, cex.lab = 1.3)
abline(h = cutoff, col = "red", lty = 2, lwd = 2)
text(obs, rhoc[obs], as.character(obs), pos = 3)

# package not longer required
detach("package:MVT")

## Fig. 2(a)
xi <- seq(0, 2, length = 200)
rho1 <- 1 - sqrt(xi)
rhoc <- 1 - xi
par(pty = "s")
plot(xi, rhoc, type = "l", xlab = "", ylab = "", xlim = c(0,2), ylim = c(-1,1), lwd = 2, lty = 2, col = "red")
lines(xi, rho1, lwd = 2)

## Fig. 2(b)
par(pty = "s")
plot(rhoc, rho1, type = "l", xlab = "", ylab = "", xlim = c(-1,1), ylim = c(-.5,1), lwd = 2)

## Fig. 3
library(plot3D)
source("bound.R")
par(pty = "s", mai = c(.5,.5,.5,.5))
persp3D(z = z, facets = FALSE, theta = 115, phi = 30, zlim = c(.5,1), xlab = "sigma11", ylab = "sigma22", zlab = "rho", bty = "b2", lwd = 2, colkey = FALSE)

# package not longer required
detach("package:plot3D")

## Table 7:

# Estimates under normal (Gaussian) distribution
fm0$center
#   manual automated 
# 2.553891  2.308982 
fm0$Scatter
#              manual automated
# manual    0.7616965 0.6942039
# automated 0.6942039 1.2369294
fm0$logLik
#[1] -200.8901

# Estimates under Laplace distribution
fm1 <- LaplaceFit(~ manual + automated, data = PSG) # converged in 18 iterations
fm1$center
#   manual automated 
# 2.569133  2.433914 
12 * fm1$Scatter
#              manual automated
# manual    0.7132149 0.6927789
# automated 0.6927789 0.9109145
fm1$logLik
#[1] -179.9035

## Table 8:

# Restricted parameter estimates under normal distribution
x <- z0$x
p <- ncol(x)
mu <- z0$Restricted$center
Sigma <- z0$Restricted$cov
u  <- chol(Sigma)
D2 <- mahalanobis(x, center = mu, cov = Sigma)
logLik.normal <- -0.5 * nobs * p * log(2 * pi) - nobs * sum(log(diag(u))) - 0.5 * sum(D2)

mu
#[1] 2.526803 2.526803
Sigma
#          [,1]      [,2]
#[1,] 0.7624303 0.6883036
#[2,] 0.6883036 1.2843752
logLik.normal
#[1] -204.7342

# Estimation of L1 coefficients (and restricted estimation)
# use larger 'bootstrap' nsamples to obtain more precision
z1 <- l1ccc(~ manual + automated, data = PSG, equal.means = TRUE, boots = TRUE, nsamples = 5000)

mu <- z1$Restricted$Fitted$center
Sigma <- z1$Restricted$Fitted$Scatter
logLik.Laplace <- z1$Restricted$Fitted$logLik

mu
#   manual automated 
# 2.529945  2.529945 
Sigma
#              manual  automated
#manual    0.06044620 0.05793709
#automated 0.05793709 0.07704239
logLik.Laplace
#[1] -183.7141

## Table 9:

# test statistics under normal distribution

center.gaussian(z0, test = "Wald")
#Wald test for the equality of location parameters 
#
#Wald statistic = 8.06, df = 1, p-value = 0.0045
#alternative hypothesis: true difference parameter is not equal to: 0 
#sample estimate: 0.2449 

center.gaussian(z0, test = "score")
#Score test for the equality of location parameters 
#
#Score statistic = 7.3387, df = 1, p-value = 0.0067
#alternative hypothesis: true difference parameter is not equal to: 0 
#sample estimate: 0.2449 

center.gaussian(z0, test = "gradient")
#Gradient test for the equality of location parameters 
#
#Gradient statistic = 7.3387, df = 1, p-value = 0.0067
#alternative hypothesis: true difference parameter is not equal to: 0 
#sample estimate: 0.2449 

# Likelihood ratio test under normal distribution
LRT0 <- -2 * (logLik.normal - fm0$logLik)
pval <- 1 - pchisq(LRT0, df = 1)
c(LRT0, pval)
#[1] 7.688086 0.005559

# test statistics under normal distribution

center.test(z1, test = "Wald")
#Wald test for the equality of location parameters 
#
#Wald statistic = 9.4267, df = 1, p-value = 0.0021
#alternative hypothesis: true difference parameter is not equal to: 0 
#sample estimate: 0.1352 
 
center.test(z1, test = "score")
#Score test for the equality of location parameters 
#
#Score statistic = 7.4012, df = 1, p-value = 0.0065
#alternative hypothesis: true difference parameter is not equal to: 0 
#sample estimate: 0.1352 
 
center.test(z1, test = "gradient")
#Gradient test for the equality of location parameters 
#
#Gradient statistic = 16.3815, df = 1, p-value = 1e-04
#alternative hypothesis: true difference parameter is not equal to: 0 
#sample estimate: 0.1352 

# Likelihood ratio test under Laplace distribution
LRT1 <- -2 * (logLik.Laplace - fm1$logLik)
pval <- 1 - pchisq(LRT1, df = 1)
c(LRT1, pval)
#[1] 7.621144 0.005769

## Table 10:

# Estimates of concordance coefficients under normal distribution
z0
#Call:
#ccc(x = ~manual + automated, data = PSG, method = "asymp", equal.means = TRUE)
#
#Coefficients:
#  estimate  variance  accuracy precision 
#   0.6744    0.0032    0.9430    0.7152 
#
#Asymptotic 95% confidence interval:
#   CCC     SE  lower  upper 
#0.6744 0.0563 0.5642 0.7847 
#
#Coefficients under equality of means:
#  estimate  variance  accuracy precision 
#   0.6726    0.0033    0.9669    0.6956 

sqrt(z0$Restricted$var.ccc) # SE of Lin's under mu_1 = mu_2
#[1] 0.057480

# L1 coefficient under mu_1 = mu_2 for normal fit
res <- unlist(z0$Restricted$L1)
res[2] <- sqrt(res[2]) # computing std.err.
res
#     rho1 var.rho1 
# 0.427780 0.050225 


# L1 estimates of concordance coefficients
z1
#Call:
#l1ccc(x = ~manual + automated, data = PSG, equal.means = TRUE, boots = TRUE, nsamples = 5000)
#
#L1 coefficients using:
#          Laplace Gaussian U-statistic
#estimate 0.5855  0.4291   0.6577     
#std.err. 0.0814  0.0937   0.0598     
#
#Lin's coefficients:
#estimate  variance  accuracy precision 
#  0.8436    0.0043    0.9815    0.8595 

# Lin's and L1 coefficients under mu_1 = mu_2 for Laplace fit
res <- unlist(z1$Restricted[c("ccc","var.ccc","rho1","var.rho1")])
res[c(2,4)] <- sqrt(res[c(2,4)]) # computing std.err.
res
#     ccc  var.ccc     rho1 var.rho1 
#0.842791 0.066949 0.603504 0.039805 
