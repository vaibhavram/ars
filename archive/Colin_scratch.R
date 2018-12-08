library(numDeriv)
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



# test with the pdf of standard normal
set.seed(2018)
f <- dnorm
D <- c(-Inf, Inf)
test_output <- get_start_points(f, D, 5)



# # param l: list of tangent lines to points in x
# # param u: list of chords between adjacent points in x
# is_reject <- function(x_star, w, l, u, z, x, h){
#   
#   j <- min(which(x_star < z))
#   u_k <- u[[j]]$intercept + u[[j]]$slope * x_star
#   i <- min(which(x_star < x))-1
#   l_k <- l[[j]]$intercept + l[[j]]$slope * x_star
#   
#   if(w <= exp(l_k - u_k)){
#     return(FALSE)
#   }
#   else{
#     if(w <= exp(h(x_star) - u_k)){
#       return(FALSE)
#     }
#     else{
#       return(TRUE)
#     }
#   }
#   
# }
# 
# is_include <- function(x_star, w, l, u, z, x, h){
#   
#   j <- min(which(x_star < z))
#   u_k <- u[[j]]$intercept + u[[j]]$slope * x_star
#   i <- min(which(x_star < x))-1
#   l_k <- l[[j]]$intercept + l[[j]]$slope * x_star
#   
#   if(w > exp(l_k - u_k)){
#     return(TRUE)
#   }
#   else{
#     return(FALSE)
#   }
#   
# }