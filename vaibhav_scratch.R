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
  return(2*x - 10*log(1 + exp(x)) - 0.5*x^2 + 100)
}

D <- c(-10, 10)

get_u_piece <- function(j, x, h, D) {
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

get_u <- function(x, h, D) {
  return(lapply(1:length(x), get_u_piece, x, h, D))
}

u <- get_u(x, h, D)

get_u_integral <- function(u, z, D) {
  full_z <- c(D[1], z, D[2])
  full_integral <- 0
  integrals <- sapply(1:length(u), function(j) {fun <- u[[j]]; (fun(full_z[j]) + fun(full_z[j+1]))/2})
  # for (j in 1:length(u)) {
  #   fun <- u[[j]]
  #   integral <- (fun(full_z[j]) + fun(full_z[j+1]))/2
  #   print(paste0(as.character(j), ": ", as.character(integral)))
  #   full_integral <- full_integral + integral
  # }
  return(integrals)
}

u_integrals <- get_u_integral(u, z, D)

# get_s <- function(u, z, D) {
#   u_integrals <- get_u_integral(u, z, D)
#   s <- lapply(1:length(u), function(t) {})
# }

sample.s <- function(n, u, z, D) {
  u_integrals <- get_u_integral(u, z, D)
  s_integrals <- u_integrals/sum(u_integrals)
  cumsum_s <- cumsum(s_integrals)
  q_vec <- runif(n)
  
  for (q in q_vec)
}




