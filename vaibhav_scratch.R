# getting the derivative of a function at a point

library(numDeriv)

# g <- function(x) {
#   return(-1*x^2 + 100)
# }
# 
# grad(g, 1:3)

#### getting u's
x <- -3:3
z <- seq(-2.5, 2.5, by = 1)

h <- function(x) {
  return(2*x - 10*log(1 + exp(x)) - 0.5*x^2 + 50)
}

D <- c(-5, 4)

# param j: index of the u piece to return
# param x: vector of x points
# param h: log of density function
# param D: domain of density function
# return: slope and intercept of tangent line at x[j]
get_u_piece <- function(j, x, h, D) {
  # throw error if j is 0 or more than k
  xj <- x[j]
  h_xj <- h(xj)
  h_prime_xj <- grad(h, xj)
  intercept <- h_xj - xj*h_prime_xj
  slope <- h_prime_xj
  return(list(intercept = intercept, slope = slope))
}

get_u <- function(x, h, D) {
  return(lapply(1:length(x), get_u_piece, x, h, D))
}

u <- get_u(x, h, D)

get_u_integral <- function(u, x, z, D) {
  full_z <- c(D[1], z, D[2])
  
  get_integral <- function(j) {
    fun <- function(t) {u[[j]]$intercept + u[[j]]$slope*t}
    return ((fun(full_z[j]) + fun(full_z[j+1]))/2)*(x[j+1] - x[j])
  }
  
  integrals <- sapply(1:length(u), get_integral)
  
  return(integrals)
}

u_integrals <- get_u_integral(u, x, z, D)

# general gist: use CDF of u to sample from s
sample.s <- function(n, u, x, z, D) {
  full_z <- c(D[1], z, D[2])
  u_integrals <- get_u_integral(u, x, z, D)
  s_integrals <- u_integrals/sum(u_integrals)
  cumsum_s <- c(0, cumsum(s_integrals))
  q <- runif(1)
  
  j <- max(which(q > cumsum_s))
  spillover <- q - cumsum_s[j]
  u_star <- u[[j]]
  z_val <- full_z[j]
  next_z <- full_z[j+1]
  
  # preparing quadratic formula
  a <- 0.5*u_star$slope
  b <- u_star$intercept
  c <- -(spillover*sum(u_integrals) + u_star$intercept*z_val + 0.5*u_star$slope*z_val^2)
  candidates <- c((-b+sqrt(b^2-4*a*c))/(2*a), (-b-sqrt(b^2-4*a*c))/(2*a))
  sample <- candidates[which(candidates < next_z & candidates > z_val)]
  return(sample)
}




