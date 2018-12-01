library(numDeriv)
library(assertthat)

#### setting up
x <- -3:3
z <- seq(-2.5, 2.5, by = 1)

h <- function(x) {
  return(2*x - 10*log(1 + exp(x)) - 0.5*x^2)
}

D <- c(-5, 4)

# param j: index of the u piece to return
# param x: vector of k points
# param h: log of density function
# return: list containing slope and intercept of tangent line at x[j]
get_u_segment <- function(j, x, h) {
  # make sure j is in range (1, k) inclusive
  assert_that(j > 0 & j <= length(x))
  
  # get h(x[j]) and h'(x[j])
  h_xj <- h(x[j])
  h_prime_xj <- grad(h, x[j])
  
  # calculate and return slope and intercept of u_segment
  intercept <- h_xj - x[j]*h_prime_xj
  slope <- h_prime_xj
  
  return(list(intercept = intercept, slope = slope))
}

# param x: vector of k points at which to find tangent lines
# param h: log of density function
# return: list with length(x) entries; jth element containing the slope
#   and intercept of the tangent line to h at x[j]
get_u <- function(x, h) {
  return(lapply(1:length(x), get_u_segment, x, h))
}

# param u: list of tangent lines to points in x
# param x: vector of k points 
# param h: log of density function
# param full_z: vector of intersection points of the tangent lines to x,
#   including domain endpoints
# return: vector of numbers, with jth element representing
#   the integral of the exponential of the tangent line of x[j] 
#   from z[j-1] to z[j] where z[0] = D[1] and z[k] = D[2]
get_s_integral <- function(u, x, h, full_z) {
  # helper function to calculate the integral of the jth
  # element of u within the proper domain
  get_integral <- function(j) {
    uj <- function(t) {u[[j]]$intercept + u[[j]]$slope*t}
    fun <- function(t) {exp(uj(t))}
    return(integrate(fun, full_z[j], full_z[j+1])$value)
  }
  
  # apply and return piecewise integrals of s
  integrals <- sapply(1:length(u), get_integral)
  return(integrals)
}

# param n: number of points to sample
# param x: vector of k points at which to find tangent lines
# param h: log of density function
# param z: vector of intersection points of the tangent lines to x
# param D: domain of density function
# return: vector of n numbers
sample.s <- function(n, x, h, z, D) {
  # combine domain endpoints with the z tangent line
  # intersection points
  full_z <- c(D[1], z, D[2])
  
  # get list of tangent lines
  u <- get_u(x, h)
  
  # get integrals under each segment of s and correspondingly of s
  s_integrals <- get_s_integral(u, x, h, full_z)
  s_integrals_norm <- s_integrals/sum(s_integrals)
  
  # create the CDF of s
  cumsum_s <- c(0, cumsum(s_integrals_norm))
  
  # draw from random uniform
  q <- runif(1) # change to n
  
  # find which u-segment q is in the domain for and
  # calculate how far into the CDF for that segment it is
  j <- max(which(q > cumsum_s))
  spillover <- q - cumsum_s[j]
  u_star <- u[[j]]
  
  # get bordering z-values for the appropriate u-segment
  z1 <- full_z[j]
  z2 <- full_z[j+1]
  
  # solve for the x* values that have the appropriate inverse CDF
  # and append to sample vector
  a <- u_star$intercept
  b <- u_star$slope
  x_star <- 1 / b *(log(b * spillover + exp(a + b * z1)) - a)
  
  return(x_star)
}


# param j: index of the l piece to return
# param x: vector of k points
# param h: log of density function
# return: list containing slope and intercept of chord from x[j] to x[j+1]
get_l_segment <- function(j, x, h) {
  # make sure j is in range (1, k - 1) inclusive
  assert_that(j > 0 & j < length(x))
  
  # solving for the common denominator
  denom <- x[j+1] - x[j]
  
  # solving for numerators for slope and intercept
  int_num <- x[j+1]*h(x[j]) - x[j]*h(x[j+1])
  slope_num <- h(x[j+1]) - h(x[j])
  
  # solving for slope and intercept
  intercept <- int_num / denom
  slope <- slope_num / denom
  return(list(intercept = intercept, slope = slope))
}

# param x: vector of k points
# param h: log of density function
# return: list with length(x) - 1 entries; jth element containing the slope
#   and intercept of the chord from x[j] to x[j+1]
get_l <- function(x, h) {
  return(lapply(1:(length(x) - 1), get_l_segment, x, h))
}






