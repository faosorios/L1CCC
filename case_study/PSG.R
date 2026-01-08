## loading packages, reading sources and attaching dataset
library(fastmatrix)
library(L1pack)
library(MVT)
data(PSG)

## computing CCC
z0 <- ccc(~ manual + automated, data = PSG, method = "asymp", equal.means = TRUE)

## removing i-th observation
nobs <- nrow(PSG)
rhoc <- rep(0, nobs)
for (i in 1:nobs)
  rhoc[i] <- ccc(~ manual + automated, data = PSG, subset = -i)$ccc

## Fig. 1
obs <- c(1,30,79)
cutoff <- z0$ccc
par(pty = "s", mai = c(1,1,.35,.35))
plot(rhoc, type = "b", ylim = c(.65, .755), ylab = "CCC estimate", lwd = 2, cex.lab = 1.3)
abline(h = cutoff, col = "red", lty = 2, lwd = 2)
text(obs, rhoc[obs], as.character(obs), pos = 3)

## model fit under normal distribution (eta = 0)
fm0 <- studentFit(~ manual + automated, data = PSG, family = Student(eta = 0))
D2 <- fm0$distances # squared Mahalanobis distances
D0 <- sqrt(D2)

## Fig. 2(a), under normality Mahalanobis distances follow a chi-distribution
cutoff <- qchi(0.975, df = 2)
par(pty = "s", mai = c(1,1,.35,.35))
plot(D0, ylim = c(0,5), ylab = "Mahalanobis distances", lwd = 2, cex.lab = 1.3)
abline(h = cutoff, lwd = 2, lty = 2, col = "red")
text(obs, D0[obs], as.character(obs), pos = 3)

## Fig. 2(b), simulated envelope under normal distribution
env <- envelope.student(fm0, reps = 5000)
ylim <- c(-2.5, 4.15)
par(pty = "s", mai = c(1,1,.35,.35))
qqnorm(env$trans, ylim = ylim, main = "", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", lwd = 2, cex.lab = 1.3)
par(new = TRUE)
qqnorm(env$envelope[,1], axes = FALSE, ylim = ylim, main = "", xlab = "", ylab = "", lwd = 2, col = "red", type = "l")
par(new = TRUE)
qqnorm(env$envelope[,2], axes = FALSE, ylim = ylim, main = "", xlab = "", ylab = "", lwd = 2, col = "red", type = "l")

## package not longer required
detach("package:MVT")

## ======================================================
## for Figures 3(a),(b) and 4, please see plots/plots.R
## ======================================================

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
z0 <- ccc(~ manual + automated, data = PSG, method = "asymp", equal.means = TRUE)
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

## reading sources for testing the equality of means
source("../code/center_test.R")

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

# test statistics under Laplace distribution

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
# ccc(x = ~manual + automated, data = PSG, method = "asymp", equal.means = TRUE)
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

## Table 11:

# removing some sets (CCC estimation)
rm.00 <- ccc(~ manual + automated, data = PSG)
rm.01 <- ccc(~ manual + automated, data = PSG, subset = -1)
rm.30 <- ccc(~ manual + automated, data = PSG, subset = -30)
rm.35 <- ccc(~ manual + automated, data = PSG, subset = -35)
rm.79 <- ccc(~ manual + automated, data = PSG, subset = -79)
rm.130 <- ccc(~ manual + automated, data = PSG, subset = -c(1,30))
rm.179 <- ccc(~ manual + automated, data = PSG, subset = -c(1,79))
rm.3079 <- ccc(~ manual + automated, data = PSG, subset = -c(30,79))
rm.all <- ccc(~ manual + automated, data = PSG, subset = -c(1,30,79))

rhoc <- c(rm.00$ccc, rm.01$ccc, rm.30$ccc, rm.35$ccc, rm.79$ccc, rm.130$ccc, rm.179$ccc, rm.3079$ccc, rm.all$ccc)
chgc <- 100 * (rhoc - rhoc[1]) / rhoc[1]

# removing some sets (L1 estimation)
rm.00 <- l1ccc(~ manual + automated, data = PSG, boots = FALSE)
rm.01 <- l1ccc(~ manual + automated, data = PSG, boots = FALSE, subset = -1)
rm.30 <- l1ccc(~ manual + automated, data = PSG, boots = FALSE, subset = -30)
rm.35 <- l1ccc(~ manual + automated, data = PSG, boots = FALSE, subset = -35)
rm.79 <- l1ccc(~ manual + automated, data = PSG, boots = FALSE, subset = -79)
rm.130 <- l1ccc(~ manual + automated, data = PSG, boots = FALSE, subset = -c(1,30))
rm.179 <- l1ccc(~ manual + automated, data = PSG, boots = FALSE, subset = -c(1,79))
rm.3079 <- l1ccc(~ manual + automated, data = PSG, boots = FALSE, subset = -c(30,79))
rm.all <- l1ccc(~ manual + automated, data = PSG, boots = FALSE, subset = -c(1,30,79))

rho1 <- c(rm.00$rho1, rm.01$rho1, rm.30$rho1, rm.35$rho1, rm.79$rho1, rm.130$rho1, rm.179$rho1, rm.3079$rho1, rm.all$rho1)
chg1 <- 100 * (rho1 - rho1[1]) / rho1[1]

rhou <- c(rm.00$ustat$rho1, rm.01$ustat$rho1, rm.30$ustat$rho1, rm.35$ustat$rho1, rm.79$ustat$rho1, rm.130$ustat$rho1, rm.179$ustat$rho1, rm.3079$ustat$rho1, rm.all$ustat$rho1)
chgu <- 100 * (rhou - rhou[1]) / rhou[1]

tab11 <- cbind(rhoc, chgc, rho1, chg1, rhou, chgu)
tab11
#      rhoc    chgc     rho1     chg1     rhou     chgu
# 0.674441  0.00000 0.585537  0.00000 0.657722  0.00000
# 0.715232  6.04813 0.617462  5.45215 0.681472  3.61082
# 0.748774 11.02143 0.627240  7.12205 0.688929  4.74464
# 0.691157  2.47854 0.607699  3.78493 0.673376  2.37999
# 0.724016  7.35059 0.615850  5.17693 0.683475  3.91545
# 0.794671 17.82664 0.660068 12.72857 0.714347  8.60924
# 0.770129 14.18775 0.649149 10.86386 0.708783  7.76324
# 0.807841 19.77942 0.659171 12.57550 0.716754  8.97522
# 0.860272 27.55342 0.693182 18.38398 0.743910 13.10399

## removing i-th observation
rho1 <- rep(0, nobs)
for (i in 1:nobs)
  rho1[i] <- l1ccc(~ manual + automated, data = PSG, boots = FALSE, subset = -i)$rho1

## Fig. 5(a)
obs <- c(1,30,35,79)
cutoff <- z1$rho1
par(pty = "s", mai = c(1,1,.35,.35))
plot(rho1, type = "b", ylim = c(.56, .64), ylab = "L1 estimate", lwd = 2, cex.lab = 1.3)
abline(h = cutoff, col = "red", lwd = 2, lty = 2)
text(obs, rho1[obs], as.character(obs), pos = 3)

## Fig. 5(b)
env1 <- envelope.Laplace(fm1, reps = 5000)
ylim <- c(-3, 4)
par(pty = "s", mai = c(1,1,.35,.35))
qqnorm(env1$trans, ylim = ylim, main = "", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", lwd = 2, cex.lab = 1.3)
par(new = TRUE)
qqnorm(env1$envelope[,1], axes = FALSE, ylim = ylim, main = "", xlab = "", ylab = "", lwd = 2, col = "red", type = "l")
par(new = TRUE)
qqnorm(env1$envelope[,2], axes = FALSE, ylim = ylim, main = "", xlab = "", ylab = "", lwd = 2, col = "red", type = "l")

## Fig. 6(a)
D1  <- fm1$distances
wts <- fm1$weights
obs <- c(1,35,79)
par(pty = "s", mai = c(1,1,.35,.35))
plot(D1, wts, xlab = "Mahalanobis distances", ylab = "Estimated weights", xlim = c(0,30), ylim = c(0,.5), lwd = 2, cex.lab = 1.3)
text(D1[obs], wts[obs], as.character(obs), pos = 3)
text(D1[30], wts[30], as.character(obs), pos = 3)

## Fig. 6(b)
obs <- c(1,30,35,79)
cutoff <- qgamma(0.976, shape = 2, scale = 2)
par(pty = "s", mai = c(1,1,.35,.35))
plot(D1, ylab = "Mahalanobis distances", ylim = c(0,30), lwd = 2, cex.lab = 1.3)
abline(h = cutoff, lwd = 2, lty = 2, col = "red")
text(obs, D1[obs], label = as.character(obs), pos = 3)

## removing i-th observation
rhou <- rep(0, nobs)
for (i in 1:nobs)
  rhou[i] <- l1ccc(~ manual + automated, data = PSG, boots = FALSE, subset = -i)$ustat$rho1

## Fig. 7
par(pty = "s", mai = c(1,1,.35,.35))
plot(rhou, type = "b", ylim = c(.645,.7), ylab = "U-stat based coefficient", lwd = 2, cex.lab = 1.3)
abline(h = rhou.00, lwd = 2, lty = 2, col = "red")
text(obs, rhou[obs], label = as.character(obs), pos = 3)
