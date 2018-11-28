# getting the derivative of a function at a point

library(numDeriv)

g <- function(x) {
  return(-1*x^2 + 100)
}

grad(g, 1:3)

#### getting u's
x <- -5:5
z <- seq(-4.5, 4.5, by = 1)

h <- function(x) {
  return(log(g(x)))
}

D <- c(-10, 10)

get_u_j <- function(j, x, z, h, D) {
  u <- function(t) {
    if (j == 0) {
      xj <- D[1]
    } else if (j == length(x)) {
      xj <- D[2]
    } else {
      xj <- x[j]
    }
    h_xj <- h(xj)
    h_prime_xj <- grad(h, xj)
    a <- h_xj - xj*h_prime_xj
    b <- h_prime_xj
    return(a+b*t)
  }
  return(u)
}

get_u <- function(x, z, h, D) {
  return(lapply(1:length(x), get_u_j, x, z, h, D))
}

u <- get_u(x, z, h, D)






