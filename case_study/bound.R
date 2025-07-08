## upper bound in Proposition 4.
bound <- function(x, y) {
  # x: sigma_11, y: sigma_22
  val <- pmin(1., 0.375 * (x + y) / sqrt(x * y))
  val
}

s1 <- seq(from = 0, to = 5, length = 26)
s1 <- s1[-1]
s2 <- seq(from = 0, to = 5, length = 26)
s2 <- s2[-1]
z <- outer(s1, s2, bound)
