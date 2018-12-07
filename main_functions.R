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
    if ( (grad(fun,x)/fun(x)) > 0 ) {
      min = x
    }
    else {
      # find the lower bound iteratively
      while( (grad(fun,x)/fun(x)) <= 0 ){
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
    if ( (grad(fun,x)/fun(x)) < 0 ) {
      max = x
    }
    else {
      # find the upper bound iteratively
      while( (grad(fun,x)/fun(x)) >= 0 ){
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
  h_prime_xj <- grad(h,x[j],method='simple')
  
  # evaluate h and h' at x[j+1]
  h_xjnext <- h(x[j+1])
  h_prime_xjnext <- grad(h, x[j+1], method='simple')
  
  # calculate numerator and denominator
  z_numerator <- h_xjnext - h_xj - ( x[j+1] * h_prime_xjnext ) + (x[j] * h_prime_xj)
  z_denominator <- h_prime_xj - h_prime_xjnext
  
  # print(z_denominator)
  
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
  assert_that(j > 0 & j <= length(x))
  
  # get h(x[j]) and h'(x[j])
  h_xj <- h(x[j])
  h_prime_xj <- grad(h, x[j])
  
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
  assert_that(j > 0 & j < length(x))
  
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

# param u: list of tangent lines to points in x
# param x: vector of k abscissae 
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
    # get limits for integral and ensure that they are different
    z1 <- full_z[j]
    z2 <- full_z[j+1]
    if (z2 < z1) {
      print(x)
      print(full_z)
      print(paste0("Z1: ", z1, " | Z2: ", z2))
    }
    
    assert_that(z2 > z1)

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
  s_integrals <- get_s_integral(u, x, h, full_z)
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
  assert_that(all(x_stars > z1s & x_stars < z2s))

  # return all x*
  return(x_stars)
}

# param j: current iteration, within D
# param x: all X for given density 
# param h: log of density function
# return: whether or not provided index satisfies upper and lower bound requirements regarding log-concave density
check_log_concave <- function(j,x,h) {
        
        # integrate u(x) over x_{index} and x_{index + 1}
        get_u_integral <- function() {
                uj <- function(t) {
                    tmp_u <- get_u_segment(j,x,h)
                    return(tmp_u$intercept + tmp_u$slope*t)
                }
                fun_u <- function(t) { exp(uj(t)) }
                return(integrate(fun_u, x[j], x[j+1])$value)
        }
        u_integral <- get_u_integral() ## store result, I_u(x)
        
        # integrate h(x) over x_{index} and x_{index + 1}
        get_h_integral <- function() {
                fun_h <- function(t) { exp(h(t)) }
                
                return(integrate(fun_h, x[j], x[j+1])$value)
        }
        h_integral <- get_h_integral() ## store result, I_h(x)
        
        # OUTPUT I_u(x) and I_h(x)
        cat("u_integral: ", u_integral, "\nh_integral: ", h_integral, "\n")
        
        # CHECK if I_u(x) > I_h(x) and OUTPUT 
        cat(ifelse(u_integral > h_integral, "upper bound segment > h(x)\n" , "ERROR upper bound segment < h(x)\n"))
        cat("\n")
        
        # ---------------------- # 
        # ---------------------- # 
        # ---------------------- # 
        
        # integrate l(x) over x_{index} and x_{index + 1}
        get_l_integral <- function() {
                lj <- function(t) {
                        tmp_l <- get_l_segment(j,x,h)
                        return(tmp_l$intercept + tmp_l$slope*t)
                }
                fun_l <- function(t) { exp(lj(t)) }
                return(integrate(fun_l, x[j], x[j+1])$value)
        }
        l_integral <- get_l_integral() ## store result, I_l(x)
        
       # OUTPUT I_l(x) and I_h(x) 
       cat("l_integral: ", l_integral, "\nh_integral: ", h_integral, "\n")
       
       # CHECK if I_l(x) < I_h(x) and OUTPUT 
       cat(ifelse(l_integral < h_integral, "lower bound segment < h(x)\n" , "ERROR lower bound segment > h(x)\n"))
       
       
       # ASSERT proper upper and lower bounds for u(x) and l(x), respectively
       # ... otherwise, error message and exit 
       assert_that(l_integral < h_integral, msg = "lower bound segment is not less than INPUT density")
       assert_that(u_integral > h_integral , msg = "upper bound segment is not greater than INPUT density")
}

# param FUN: density function from which to sample
# param n: number of points to sample
# param D: domain of density function, a numeric vector of
#   length two
# param verbose (optional): verbose output desired?
# return: n points sampled from FUN using adaptive-rejection 
#   sampling
ars <- function(FUN, n = 1, D = c(-Inf, Inf), verbose = FALSE){
  
  # checking classes for each argument
  assert_that(class(n) == "numeric")
  assert_that(class(FUN) == "function")
  assert_that(class(D) == "numeric")
  
  # assure that n is admissible
  assert_that(length(n) == 1, n > 0)
  
  # assure that D is admissible
  assert_that(length(D) == 2, D[2] > D[1])
  
  # normalize function
  fun_integral <- integrate(FUN, D[1], D[2])
  assert_that(fun_integral$value > 0)
  f <- function(t) FUN(t)/fun_integral$value
  
  # initialize sample
  sample <- c()
  
  # get h(t) = log(f(t))
  h <- function(x) {
    return(log(f(x)))
  }
  
  # initialize abscissae and batch.size
  x <- get_start_points(f, D)
  #x <- c(1, 2, 3)
  batch.size <- 1
  
  # index to track iterations, if verbose=TRUE
  i <- 1
  
  while(length(sample) < n){
    
    # print out info if verbose
    if (verbose) {
      cat("Batch ", i, ":\n", sep = "")
    }
    
    # update z and make sure z corresponds in dimension
    z <- get_z_all(x, h)
    assert_that(length(z) + 1 == length(x))
    
    # combine domain endpoints with the z tangent line
    # intersection points
    full_z <- c(D[1], z, D[2])
    
    # print(full_z)
    
    # get upper bound and lower bound
    u <- get_u(x, h)
    l <- get_l(x, h)
    
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
    x <- sort(c(x, x_stars[!(check1)]))
    
    # print info if verbose
    if (verbose) {
      cat(" Batch Size:", batch.size, "\n")
      cat(" Accepted:", sum(check1 | check2), "\n")
      cat(" Rejected:", batch.size - sum(check1 | check2), "\n")
      cat(" Failed Check 1:", sum(!check1), "\n")
      i <- i + 1
    }
    
    # increases batch size if none failed check 1
    if (all(check1)) {
      batch.size <- 2 * batch.size
    }
  }
  
  # return sample of length n
  return(sample[1:n])
}


