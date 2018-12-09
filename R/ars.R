library(numDeriv)
library(assertthat)

# param fun: the density function
# param D: domain of density function
# param n: number of starting points
# return: list containing n starting points
get_start_points <- function(fun, D, n=3, x_start=2, x_step=1){
  
  # check if the lower bound is finite
  if (is.finite(D[1])){
    min = D[1] + 0.001
  }
  else{
    x = -x_start
    # if its first derivative is negative, keep it as the lower bound
    if ( (numDeriv::grad(fun,x)/fun(x)) > 0 ) {
      min = x
    }
    else {
      # find the lower bound iteratively
      while( (numDeriv::grad(fun,x)/fun(x)) <= 0 ){
        x <- x-x_step
      }
      min = x 
    }
  }
  
  # check if the upper bound is finite  
  if (is.finite(D[2])){
    max = D[2] - 0.001
  }
  else {
    x = x_start
    # if its first derivative is negative, keep it as the upper bound
    if ( (numDeriv::grad(fun,x)/fun(x)) < 0 ) {
      max = x
    }
    else {
      # find the upper bound iteratively
      while( (numDeriv::grad(fun,x)/fun(x)) >= 0 ){
        x <- x+x_step
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


# param j: index of abscissae for which to find
#   the tangent line intersection point
# param x: vector of k abscissae
# param h: log of density function
# return: intersection of lines tangent to h
#   at x[j] and x[j+1]
get_z <- function(j, x, h, eps = 1e-08) {
  
  # evaluate h and h' at x[j]
  h_xj <- h(x[j])
  h_prime_xj <- numDeriv::grad(h,x[j],method='simple')
  
  # evaluate h and h' at x[j+1]
  h_xjnext <- h(x[j+1])
  h_prime_xjnext <- numDeriv::grad(h, x[j+1], method='simple')
  
  # calculate numerator and denominator
  z_numerator <- h_xjnext - h_xj - ( x[j+1] * h_prime_xjnext ) + (x[j] * h_prime_xj)
  z_denominator <- h_prime_xj - h_prime_xjnext
  
  # if h is a straight line, just get average of x's
  if (abs(z_denominator) < eps) {
    return((x[j]+x[j+1])/2)
  } else {
    return(z_numerator / z_denominator)  
  }
}


# param x: vector of k abscissae
# param h: log of density function
# param D: domain of density function
# return: vector of k - 1 tangent line
#   intersection points
get_z_all <- function(x, h) {
  return(sapply(1:(length(x) - 1), get_z, x, h))
}


# param j: index of the u piece to return
# param x: vector of k abscissae
# param h: log of density function
# return: list containing slope and intercept of tangent line at x[j]
get_u_segment <- function(j, x, h) {
  # make sure j is in range (1, k) inclusive
  assertthat::assert_that(j > 0 & j <= length(x))
  
  # get h(x[j]) and h'(x[j])
  h_xj <- h(x[j])
  h_prime_xj <- numDeriv::grad(h, x[j])
  
  # calculate and return slope and intercept of u_segment
  return(list(intercept = h_xj - x[j]*h_prime_xj, slope = h_prime_xj))
}

# param x: vector of k abscissae
# param h: log of density function
# return: list of k entries; jth element containing the slope
#   and intercept of the tangent line to h at x[j]
get_u <- function(x, h) {
  return(lapply(1:length(x), get_u_segment, x, h))
}



# param j: index of the l segment to return
# param x: vector of k abscissae
# param h: log of density function
# return: list containing slope and intercept of chord 
#   from x[j] to x[j+1]
get_l_segment <- function(j, x, h) {
  # make sure j is in range (1, k - 1) inclusive
  assertthat::assert_that(j > 0 & j < length(x))
  
  # solving for the common denominator
  denom <- x[j+1] - x[j]
  
  # solving for numerators for slope and intercept
  int_num <- x[j+1]*h(x[j]) - x[j]*h(x[j+1])
  slope_num <- h(x[j+1]) - h(x[j])
  
  # calculate and return slope and intercept
  return(list(intercept = int_num / denom, slope = slope_num / denom))
}

# param x: vector of k abscissae
# param h: log of density function
# return: list with k - 1 entries; jth element containing the slope
#   and intercept of the chord from x[j] to x[j+1]
get_l <- function(x, h) {
  return(lapply(1:(length(x) - 1), get_l_segment, x, h))
}

# param l: lower hull; list of chords where jth element is slope
#   and intercept of line from x[j] to x[j+1]
# param x: vector from original k abscissae
# return: vector of l integrals , for comparison with density and check
#    for log concavity
get_l_integral <- function(l, x) {
  
  get_integral <- function(j) {
    
    # store initial x 
    #    and next x
    x_first <- x[j]
    x_next <- x[j+1]
    
    # store intercept and slope at initial x
    a <- l[[j]]$intercept
    b <- l[[j]]$slope
    
    # integral (analytical, via Vaibhav)
    if (b == 0) {
      return(exp(a) * (x_next - x_first))
    } else {
      return(1 / b * (exp(a + b * x_next) - exp(a + b * x_first)))
    }
  }
  
  # iterate through all of relevant x in l)
  l_integral_all <- sapply(1:length(l), get_integral)
  
  # return vector of integrals
  #    for downstream comparison with density 
  return(l_integral_all)
}

# param u: list of tangent lines to points in x
# param full_z: vector of intersection points of the tangent lines to x,
#   including domain endpoints
# return: vector of numbers, with jth element representing
#   the integral of the exponential of the tangent line of x[j] 
#   from z[j-1] to z[j] where z[0] = D[1] and z[k] = D[2]
get_s_integral <- function(u, full_z) {
  # helper function to calculate the integral of the jth
  # element of u within the proper domain
  get_integral <- function(j) {
    # get limits for integral and ensure that z2 > z1
    z1 <- full_z[j]
    z2 <- full_z[j+1]
    
    assertthat::assert_that(z2 > z1)
    
    # get slope and intercept of u[[j]]
    a <- u[[j]]$intercept
    b <- u[[j]]$slope
    
    # return the analytical solution to the integral
    if (b == 0) {
      return(exp(a) * (z2 - z1))
    } else {
      return(1 / b * (exp(a + b * z2) - exp(a + b * z1)))
    }
  }
  
  # apply and return piecewise integrals of s
  integrals <- sapply(1:length(u), get_integral)
  return(integrals)
}

# param f: a density function
# param vec: a vector of points of length k that define intervals
#   over which to integrate
# return: a vector of integrals of length k-1 in which the jth
#   element is the integral under f between vec[j] and vec[j+1]
get_f_integral <- function(f, vec) {
  # helper function to calculate the integral under f from
  # vec[j] to vec[j+1]
  get_integral <- function(j) {
    # get limits for integral and ensure that z2 > z1
    a <- vec[j]
    b <- vec[j+1]
    assertthat::assert_that(b > a)
    
    # return the integral
    return(integrate(f, a, b)$value)
  }
  
  # apply and return piecewise integrals of s
  integrals <- sapply(1:(length(vec) - 1), get_integral)
  return(integrals)
}

# param n: number of points to sample
# param x: vector of k abscissae
# param h: log of density function
# param full_z: vector of intersection points of the tangent 
#   lines to x, with domain endpoints as well
# param u: piecewise upper bound fn for h
# return: vector of n numbers sampled from S,
#   where S = u / integral(u)
sample.s <- function(n, x, h, full_z, u) {
  
  # get integrals under each segment of s
  # and full integral under s
  s_integrals <- get_s_integral(u, full_z)
  full_s_integral <- sum(s_integrals)
  
  # get normalized integrals under s
  s_integrals_norm <- s_integrals/full_s_integral
  
  # create the CDF of s
  cumsum_s <- c(0, cumsum(s_integrals_norm))
  
  # draw from random uniform
  qs <- runif(n)
  
  # get index of u_segments for which each q is in
  # and calculate how far into the CDF for that
  # segment it is
  js <- sapply(qs, function(q) max(which(q > cumsum_s)))
  spillovers <- qs - cumsum_s[js]
  
  # get intervals for the domain in which each x* should fall
  z1s <- full_z[js]
  z2s <- full_z[js+1]
  
  # get  u* for each x* and corresponding slope and intercept
  u_stars <- lapply(js, function(j) u[[j]])
  as <- sapply(1:n, function(i) u_stars[[i]]$intercept)
  bs <- sapply(1:n, function(i) u_stars[[i]]$slope)
  
  # get each x*, using analytical solution to integral
  x_stars <- sapply(1:n, function(i) {
    if(bs[i] != 0){
      return((log(bs[i] * spillovers[i] * full_s_integral + exp(as[i] + bs[i] * z1s[i])) - as[i]) / bs[i])
    }
    else{
      return((spillovers[i] * full_s_integral)/exp(as[i]) + z1s[i])
    }
  })
  
  # ensure that all x*s are in proper domain
  assertthat::assert_that(all(x_stars > z1s & x_stars < z2s))
  
  # return all x*
  return(x_stars)
}

# param fun: a log density
# param D: the domain of the density
# param eps (optional): threshhold value for difference
#     comparisons
# return: TRUE if fun is a linear function
is_linear <- function(fun, D, eps = 1e-08){
  
  # check if domain limits are finite and if not,
  # set manually
  if(is.finite(D[1])){
    min <- D[1]
  } else {
    min <- min(-100, D[2]-1)
  }
  
  if(is.finite(D[2])){
    max <- D[2]
  } else {
    max <- max(100, D[1]+1)
  }
  
  # sample 100 test points to apply gradient
  # and keep only those for which fun(t) is finite
  test <- runif(100, min, max)
  test <- test[is.finite(fun(test))]
  
  # apply gradient to elements of test
  # and check and return if they are all the same
  results <- numDeriv::grad(fun, test)
  results <- results[!is.na(results)]
  differences <- abs(results - mean(results)) < eps
  return(all(differences))
}

# param FUN: density function from which to sample
# param n: number of points to sample
# param D: domain of density function, a numeric vector of
#   length two
# return: n points sampled from FUN using adaptive-rejection 
#   sampling
ars <- function(FUN, n = 1, D = c(-Inf, Inf)){
  
  # checking classes for each argument
  assertthat::assert_that(class(n) == "numeric")
  assertthat::assert_that(class(FUN) == "function")
  assertthat::assert_that(class(D) == "numeric")
  
  # assure that n is admissible
  assertthat::assert_that(length(n) == 1, n > 0)
  
  # assure that D is admissible
  assertthat::assert_that(length(D) == 2, D[2] > D[1])
  
  # normalize function
  fun_integral <- integrate(FUN, D[1], D[2])
  assertthat::assert_that(fun_integral$value > 0)
  f <- function(t) FUN(t)/fun_integral$value
  
  # initialize sample
  sample <- c()
  
  # get h(t) = log(f(t))
  h <- function(x) {
    return(log(f(x)))
  }
  
  # is_linear has a fail rate of < 0.5% so we run it twice
  # to decrease the probability that it gives a false negative
  is_linear <- is_linear(h,D) | is_linear(h,D)
  
  # initialize abscissae and batch.size
  x <- get_start_points(f, D)
  batch.size <- 1
  
  while(length(sample) < n){
    
    # update z and make sure z corresponds in dimension
    z <- get_z_all(x, h)
    assertthat::assert_that(length(z) + 1 == length(x))
    
    # combine domain endpoints with the z tangent line
    # intersection points
    full_z <- c(D[1], z, D[2])
    
    # get upper bound and lower bound
    u <- get_u(x, h)
    l <- get_l(x, h)
    
    # don't need to check concavity if it's linear
    if(! is_linear){
      # running concavity check for upper bound
      s_integrals <- get_s_integral(u, full_z)
      f_integrals_z <- get_f_integral(f, full_z)
      assertthat::assert_that(all(s_integrals >= f_integrals_z))
      
      # running concavity check for lower bound
      l_integrals <- get_l_integral(l, x)
      f_integrals_x <- get_f_integral(f, x)
      assertthat::assert_that(all(l_integrals <= f_integrals_x))
    }
    
    # get sample of size batch.size from s
    x_stars <- sample.s(batch.size, x, h, full_z, u)
    
    # draw sample from Unif(0,1) of size batch.size
    ws <- runif(batch.size)
    
    # get index of corresponding u segment for each x*
    js <- sapply(x_stars, function(x_star) min(which(x_star < full_z)))
    
    # evaluate U_k and L_k at each x*
    uks_xstar <- sapply(1:batch.size, function(i) u[[js[i]-1]]$intercept + u[[js[i]-1]]$slope * x_stars[i])
    lks_xstar <- sapply(1:batch.size, function(index) {
      if (x_stars[index] > x[length(x)] || x_stars[index] < x[1]) {
        return(-Inf)
      } else {
        i <- min(which(x_stars[index] < x))-1
        return(l[[i]]$intercept + l[[i]]$slope * x_stars[index])
      }
    })
    
    # check if w <= exp(L(x*) - U(x*)) for each x*
    check1 <- ws <= exp(lks_xstar - uks_xstar)
    
    # check if w <= exp(h(x*) - U(x*)) for each x*
    check2 <- ws <= exp(h(x_stars) - uks_xstar)
    
    # add points x* which pass check 1 or 2 to the sample
    sample <- c(sample, x_stars[check1 | check2])
    
    # add points x* which fail check 1 to vector of abscissae
    # and sort
    if(! is_linear){
      x <- sort(c(x, x_stars[!(check1)]))
    }
    
    # increases batch size if none failed check 1
    if (all(check1)) {
      batch.size <- 2 * batch.size
    }
  }
  
  # return sample of length n
  return(sample[1:n])
}

