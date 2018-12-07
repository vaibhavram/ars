library(testthat)

source("main_functions.R")

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

test_that("Run different tests and compare the results",{
  significance<-0.05
  x_norm<-ars(dnorm, n=10000, D=c(-Inf, Inf), verbose = FALSE)
  output_norm<-ks.test(x_norm, pnorm)
  expect_lt(significance, output_norm$p.value)
  x_exp<-ars(dexp, n=10000, D=c(0, Inf), verbose = FALSE)
  output_exp<-ks.test(x_exp, pexp)
  expect_lt(significance, output_exp$p.value)
  x_beta<-ars(function(t) dbeta(t, 2, 2), n=10000, D=c(0, 1), verbose = FALSE)
  output_beta<-ks.test(x_beta, pbeta, 2, 2)
  expect_lt(significance, output_beta$p.value)
  x_gamma<-ars(function(t) dgamma(t, 2), n=10000, D=c(0, Inf), verbose = FALSE)
  output_gamma<-ks.test(x_gamma, pgamma, 2)
  expect_lt(significance, output_gamma$p.value)
  x_chisq2<-ars(function(t) dchisq(t, 2), n=10000, D=c(0, Inf), verbose = FALSE)
  output_chisq2<-ks.test(x_chisq2, pchisq, 2)
  expect_lt(significance, output_chisq2$p.value)
  x_chisq3<-ars(function(t) dchisq(t, 3), n=10000, D=c(0, Inf), verbose = FALSE)
  output_chisq3<-ks.test(x_chisq3, pchisq, 3)
  expect_lt(significance, output_chisq3$p.value)
  x_unif<-ars(dunif, n=10000, D=c(0, 1), verbose = FALSE)
  output_unif<-ks.test(x_unif, punif)
  expect_lt(significance, output_unif$p.value)
  })

