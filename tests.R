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

test_that("get_u() returns corrent slopes of tangent lines for simple h", {
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




