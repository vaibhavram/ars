library(numDeriv)
library(assertthat)
library(testthat)

# param fun: the density function
# param n: number of starting points
# param D: domain of density function
# return: list containing n starting points
# we assume that the input function is log concave
get_start_points <- function(fun, D, n=3){
  
  # check if the lower bound is finite
  if (is.finite(D[1])){
    min = D[1]
  }
  else{
    x = -10
    # if its first derivative is negative, keep it as the lower bound
    if ( (grad(fun,x)/fun(x)) > 0 ) {
      min = x
    }
    else {
      # find the lower bound iteratively
      while( (grad(fun,x)/fun(x)) <= 0 ){
        x <- x-1
      }
      min = x 
    }
  }
  
  # check if the upper bound is finite  
  if (is.finite(D[2])){
    max = D[2]
  }
  else {
    x = 10
    # if its first derivative is negative, keep it as the upper bound
    if ( (grad(fun,x)/fun(x)) < 0 ) {
      max = x
    }
    else {
      # find the upper bound iteratively
      while( (grad(fun,x)/fun(x)) >= 0 ){
        x <- x+1
      }
      max =x
    }
  }
  
  m = n-2
  # generate the rest numbers from uniform dist
  nums <- runif(m, min, max)
  output <- c(min,sort(nums),max)
  
  return(output)
}


#### FROM BRANDON - CHANGE IF HE CHANGES

get_z  <- function(x_initial, x_next, h) {
  h_xj <- h(x_initial)
  h_prime_xj <- grad(h,x_initial)
  
  h_xjnext <- h(x_next)
  h_prime_xjnext <- grad(h, x_next)
  
  z_numerator <- h_xjnext - h_xj - ( x_next * h_prime_xjnext ) + (x_initial * h_prime_xj)
  z_denominator <- h_prime_xj - h_prime_xjnext
  return( z_numerator / z_denominator )
}

get_z_all  <- function(x, h, D) {
  k <- length(x)
  store_all_z <- c()
  for (i in 1:(k-1)) {
    store_all_z <- c(store_all_z, get_z(x[i], x[i+1], h))
  }
  return(store_all_z)
}



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
sample.s <- function(x, h, z, D) {
  # check to make sure h is concave
  #assert_that(is.concave(h)) # requires function that brandon will write
  
  # make sure we are sampling at least 1 point
  #assert_that(n > 0)
  
  # make sure z and x correspond in dimension
  assert_that(length(z) + 1 == length(x))
  
  # make sure D is properly formatted
  assert_that(length(D) == 2)
  
  # combine domain endpoints with the z tangent line
  # intersection points
  full_z <- c(D[1], z, D[2])
  
  # get list of tangent lines
  u <- get_u(x, h)
  
  # get integrals under each segment of s
  s_integrals <- get_s_integral(u, x, h, full_z)
  denom <- sum(s_integrals)
  # make sure the total integral under s is more than the 
  # integral under g (since s is an upper bound)
  test_int <- integrate(function(t) exp(h(t)), D[1], D[2])
  assert_that(denom > test_int$value - test_int$abs.error)
  # get normaized integrals under s
  s_integrals_norm <- s_integrals/denom
  
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
  x_star <- (log(b * spillover + exp(a + b * z1)) - a) / b
  
  # make sure x* is between z1 and z2
  assert_that(x_star >= z1, x_star <= z2)
  
  return(x_star)
}

# param l: list of tangent lines to points in x
# param u: list of chords between adjacent points in x
is_reject <- function(x_star, w, l, u, z, x, h){
  
  j <- min(which(x_star < z))
  u_k <- u[[j]]$intercept + u[[j]]$slope * x_star
  i <- min(which(x_star < x))-1
  l_k <- l[[j]]$intercept + l[[j]]$slope * x_star
  
  if(w <= exp(l_k - u_k)){
    return(FALSE)
  }
  else{
    if(w <= exp(h(x_star) - u_k)){
      return(FALSE)
    }
    else{
      return(TRUE)
    }
  }
  
}

is_include <- function(x_star, w, l, u, z, x, h){
  
  j <- min(which(x_star < z))
  u_k <- u[[j]]$intercept + u[[j]]$slope * x_star
  i <- min(which(x_star < x))-1
  l_k <- l[[j]]$intercept + l[[j]]$slope * x_star
  
  if(w > exp(l_k - u_k)){
    return(TRUE)
  }
  else{
    return(FALSE)
  }
  
}
  

#### setting up
set.seed(736)
D <- c(-Inf, Inf)


### sanity check
# x <- get_start_points(fun, D, n=5)
# z <- get_z_all(x, h, D)
# u <- get_u(x, h)
# l <- get_l(x, h)
# x_star <- sample.s(x, h, z, D)
# is_reject(x_star, l, u)


ars <- function(n, fun, D){
  sample = c()
  h <- function(x) {
    return(log(fun(x)))
  }
  x <- get_start_points(f, D)
  
  while(length(sample) < n){
    z <- get_z_all(x, h, D)
    u <- get_u(x, h)
    l <- get_l(x, h)
    x_star <- sample.s(x, h, z, D)
    
    # sampling step
    w <- runif(1)
    j <- min(which(x_star < z))
    u_k <- u[[j]]$intercept + u[[j]]$slope * x_star
    i <- min(which(x_star < x))-1
    l_k <- l[[j]]$intercept + l[[j]]$slope * x_star
    
    if(w <= exp(l_k - u_k)){
      sample = c(sample, x_star)
    }
    else{
      # Updating step
      x <- sort(c(x, x_star))
      if(w <= exp(h(x_star) - u_k)){
        sample = c(sample, x_star)
      }
    }
  }
  return(sample)
}
