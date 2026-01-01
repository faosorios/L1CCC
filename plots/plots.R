## required package for persp3D function
library(plot3D)

## Fig. 3(a)
xi <- seq(0, 2, length = 200)
rho1 <- 1 - sqrt(xi)
rhoc <- 1 - xi
par(pty = "s")
plot(xi, rhoc, type = "l", xlab = expression(xi), ylab = "agreement coefficient", xlim = c(0,2), ylim = c(-1,1), lwd = 2, lty = 2, col = "red", cex.lab = 1.3)
lines(xi, rho1, lwd = 2)

## Fig. 3(b)
par(pty = "s")
plot(rhoc, rho1, type = "l", xlab = expression(rho[c]), ylab = expression(rho[1]), xlim = c(-1,1), ylim = c(-.5,1), lwd = 2, cex.lab = 1.3)

## Fig. 4
source("bound.R")
par(pty = "s", mai = c(.5,.5,.5,.5))
persp3D(z = z, facets = FALSE, theta = 115, phi = 30, zlim = c(.5,1), xlab = "sigma11", ylab = "sigma22", zlab = "rho", bty = "b2", lwd = 2, colkey = FALSE)

# package not longer required
detach("package:plot3D")