

### VAIBHAV

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