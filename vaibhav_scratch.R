# getting the derivative of a function at a point

library(numDeriv)

g <- function(x) {
  return(-2 * x^2 + 100)
}

grad(g, 1:3)


#### getting u's
x <- -5:5
z <- seq(-4.5, 4.5, by = 1)

h <- function(x) {
  return(log(g(x)))
}

get_u_j <- function(j, x, z, h) {
  
}









