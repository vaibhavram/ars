library(testthat)

source("ars.R")

context("Tests for ars package")

test_that("get_u_segment() returns correct tangent line for simple h", {
  j <- 1
  x <- c(0)
  offset <- 5
  h <- function(t) -t^2 + offset
  u_seg <- get_u_segment(j, x, h)
  expect_equal(u_seg$slope, 0)
  expect_equal(u_seg$intercept, offset)
})

test_that("get_u() returns correct slopes of tangent lines for simple h", {
  x <- c(-1, 1)
  h <- function(t) -t^2
  u <- get_u(x, h)
  u1 <- u[[1]]
  u2 <- u[[2]]
  expect_equal(u1$slope, 2)
  expect_equal(u2$slope, -2)
})

test_that("get_s_integral() returns appropriate upper-bound integral", {
  x <- seq(-4.95, 4.95, length.out = 200)
  h <- function(t) -t^2
  u <- get_u(x, h)
  D <- c(-5, 5)
  z <- sapply(1:(length(x) - 1), function(i, vec) (vec[i] + vec[i+1])/2, x)
  full_z <- c(D[1], z, D[2])
  s_ints <- get_s_integral(u, z, h, full_z)
  s_int <- sum(s_ints)
  tru_int <- integrate(f = function(t) exp(h(t)), D[1], D[2])
  epsilon <- 0.01
  expect_lt(s_int - tru_int$value, epsilon)
})

test_that("Sampling Normal(0,1) through ars() passes K-S test",{
  significance <- 0.05
  
  x_norm <- ars(dnorm, n=10000, D=c(-Inf, Inf), verbose = FALSE)
  output_norm <- ks.test(x_norm, pnorm)
  expect_lt(significance, output_norm$p.value)
})

test_that("Sampling Normal(1000,1) through ars() passes K-S test",{
  significance <- 0.05
  
  x_norm <- ars(function(t) dnorm(t, 1000, 1), n=10000, D=c(-Inf, Inf), verbose = FALSE)
  output_norm <- ks.test(x_norm, pnorm, 1000, 1)
  expect_lt(significance, output_norm$p.value)
})

test_that("Sampling Beta(1,1) through ars() passes K-S test",{
  significance <- 0.05
  
  x_beta <- ars(function(t) dbeta(t, 1, 1), n=10000, D=c(0, 1), verbose = FALSE)
  output_beta <- ks.test(x_beta, pbeta, 1, 1)
  expect_lt(significance, output_beta$p.value)
})

test_that("Sampling Beta(2,2) through ars() passes K-S test",{
  significance <- 0.05
  
  x_beta2 <- ars(function(t) dbeta(t, 2, 2), n=10000, D=c(0, 1), verbose = FALSE)
  output_beta2 <- ks.test(x_beta2, pbeta, 2, 2)
  expect_lt(significance, output_beta2$p.value)
})

test_that("Sampling Gamma(1, 1) through ars() passes K-S test",{
  significance <- 0.05
  
  x_gamma <- ars(function(t) dgamma(t, 1), n=10000, D=c(0, Inf), verbose = FALSE)
  output_gamma <- ks.test(x_gamma, pgamma, 1)
  expect_lt(significance, output_gamma$p.value)
})

test_that("Sampling Gamma(2, 1) through ars() passes K-S test",{
  significance <- 0.05
  
  x_gamma2 <- ars(function(t) dgamma(t, 2), n=10000, D=c(0, Inf), verbose = FALSE)
  output_gamma2 <- ks.test(x_gamma2, pgamma, 2)
  expect_lt(significance, output_gamma2$p.value)
})

test_that("Sampling ChiSq(2) through ars() passes K-S test",{
  significance <- 0.05
  
  x_chisq2 <- ars(function(t) dchisq(t, 2), n=10000, D=c(0, Inf), verbose = FALSE)
  output_chisq2 <- ks.test(x_chisq2, pchisq, 2)
  expect_lt(significance, output_chisq2$p.value)
})

test_that("Sampling ChiSq(3) through ars() passes K-S test",{
  significance <- 0.05
  
  x_chisq3 <- ars(function(t) dchisq(t, 3), n=10000, D=c(0, Inf), verbose = FALSE)
  output_chisq3 <- ks.test(x_chisq3, pchisq, 3)
  expect_lt(significance, output_chisq3$p.value)
})

test_that("Sampling Exp(1) through ars() passes K-S test",{
  significance <- 0.05
  
  x_exp <- ars(dexp, n=10000, D=c(0, Inf), verbose = FALSE)
  output_exp <- ks.test(x_exp, pexp)
  expect_lt(significance, output_exp$p.value)
})

test_that("Sampling Unif(0,1) through ars() passes K-S test",{
  significance <- 0.05
  
  x_unif <- ars(dunif, n=10000, D=c(0, 1), verbose = FALSE)
  output_unif <- ks.test(x_unif, punif)
  expect_lt(significance, output_unif$p.value)
})

test_that("Sampling t(df = 3) throws an error, since it is not log-concave", {
  f <- function(t) dt(t, df = 3)
  expect_error(ars(f, 1000))
  expect_error(ars(f, 1))
})

test_that("Sampling Gamma(0.5, 1) throws an error, since it is not log-concave", {
  f <- function(t) dgamma(t, 0.5)
  expect_error(ars(f, 1000))
  expect_error(ars(f, 1))
})