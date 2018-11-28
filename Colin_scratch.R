library(numDeriv)

# we assume that the input function is log concave
random_num_generator <- function(fun, k=3){

  # randomly pick a number between 0 and 1
  x <- runif(1)
  # if its first derivative is negative
  if ( (grad(fun,x)/fun(x)) < 0 ) {
    max = x
    # find the upper bound iteratively
    while( (grad(fun,x)/fun(x)) <= 0 ){
      x <- runif(1, x-1, x)
    }
    min = x
  }
  # if its first derivative is zero or positive
  else{
    if( (grad(fun,x)/fun(x)) == 0 ){
    x = x+1
    }
    min = x
    # find the upper bound iteratively
    while( (grad(fun,x)/fun(x)) >= 0){
      x <- runif(1, x, x+1)
    }
    max = x
  }

  m = k-2
  nums <- runif(m, min, max)
  output <- c(min,sort(nums),max)
  
  return(output)
}

# test with the pdf of exponential
set.seed(0)
f <- dnorm
test_output <- random_num_generator(f, k=10)