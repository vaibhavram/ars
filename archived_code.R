# lks <- sapply(x_stars, function(x_stars, x) {
#   if (x_star > x[length(x)] || x_star < x[1]){
#     l_k = -Inf
#   }
#   else{
#     i <- min(which(x_star < x))-1
#     l_k <- l[[i]]$intercept + l[[i]]$slope * x_star
#   }
# }, x)
# for (x_star in x_stars) {
#   # sampling step
#   w <- runif(1)
#   #print(paste0('w is',w))
#   #print(paste0('z is',z))
#   j <- min(which(x_star < full_z))
#   u_k <- u[[j-1]]$intercept + u[[j-1]]$slope * x_star
#   #print(paste0('U-k is',u_k))
#   if (x_star > x[length(x)] || x_star < x[1]){
#     l_k = -Inf
#   }
#   else{
#     i <- min(which(x_star < x))-1
#     l_k <- l[[i]]$intercept + l[[i]]$slope * x_star
#   }
#   #print(paste0('L-k is',l_k))
#   
#   if(w <= exp(l_k - u_k)){
#     sample = c(sample, x_star)
#   }
#   else{
#     # Updating step
#     x <- sort(c(x, x_star))
#     if(w <= exp(h(x_star) - u_k)){
#       sample = c(sample, x_star)
#     }
#   }
# }

# # get js and spillovers
# js <- sapply(qs, function(q) max(which(q > cumsum_s)))
# spillovers <- qs - cumsum_s[js]
# u_stars <- lapply(js, function(j) u[[j]])
# z1s <- full_z[js]
# z2s <- full_z[js+1]
# as <- sapply(1:batch.size, function(i) u_star[[i]]$intercept)
# bs <- sapply(1:batch.size, function(i) u_star[[i]]$slope)
# x_stars <- sapply(1:batch.size, function(i) {
#   if (bs[i] != 0) {
#     return((log(bs[i] * s_integrals_norm * spillovers[i] + exp(as[i] + bs[i] * z1s[i])) - as[i]) / bs[i])
#   }
#   else {
#     return((s_integrals_norm * spillovers[i])/exp(as[i]) + z1s[i])
#   }
# })
# assert_that(all(x_stars < z2s, x_stars > z1s))
# 
# for (q in qs) {
#   # find which u-segment q is in the domain for and
#   # calculate how far into the CDF for that segment it is
#   j <- max(which(q > cumsum_s))
#   spillover <- q - cumsum_s[j]
#   u_star <- u[[j]]
#   
#   # get bordering z-values for the appropriate u-segment
#   z1 <- full_z[j]
#   z2 <- full_z[j+1]
#   
#   # solve for the x* values that have the appropriate inverse CDF
#   # and append to sample vector
#   a <- u_star$intercept
#   b <- u_star$slope
#   if(b != 0){
#     x_star <- (log(b * s_integrals_norm * spillover + exp(a + b * z1)) - a) / b
#   }
#   else{
#     x_star <- (s_integrals_norm * spillover)/exp(a) + z1
#   }
#   
#   # make sure x* is between z1 and z2
#   # if (! assert_that(x_star >= z1, x_star <= z2)) {
#   #   print(paste0("x*: ", x_star))
#   #   print(paste0("Z: ", c(z1, z2)))
#   # }
#   
#   x_stars <- c(x_stars, x_star)
# }