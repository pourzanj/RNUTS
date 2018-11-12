library(mockery)

context("The Hamiltonian System for a Multivariate Gaussian")

######################
# Simple 2x2
######################
M <- diag(c(1.0, 1.0))
Sigma <- matrix(c(1.0, 0.9, 0.9, 1.0), nrow = 2)
gaussian_system <- create_gaussian_hamiltonian_system(M, Sigma)

z1 <- list(q = c(0.0, 0.0), p = c(-0.1759228, 0.8358274))
z2 <- list(q = c(0.08784589, -0.96785614), p = c(-0.1759228, 0.8358274))

test_that("M is included and matches input", {
  expect_equal(gaussian_system$M, M)
})

test_that("Hamiltonian is computed correctly", {
  expect_equal(gaussian_system$compute_H(z1), 0.5*sum(z1$p*solve(M,z1$p)))
  expect_equal(gaussian_system$compute_H(z2), 0.5*sum(z2$p*solve(M,z2$p)) + 0.5*sum(z2$q*solve(Sigma,z2$q)))
})

test_that("gradU is computed correctly", {
  expect_equal(gaussian_system$compute_gradU(z1$q), c(0.0, 0.0))
  expect_equal(gaussian_system$compute_gradU(z2$q), solve(Sigma,c(0.08784589, -0.96785614)))
})

test_that("hessU vector-product is computed correctly", {
  expect_equal(gaussian_system$compute_hessU_vec_prod(z1$q, c(0.0, 0.0)), c(0.0, 0.0))
  expect_equal(gaussian_system$compute_hessU_vec_prod(z2$q, c(0.5, 0.123)), solve(Sigma,c(0.5, 0.123)))
})

test_that("momentum sampler is included and draws from the correct momentum", {
  L <- chol(M)
  mock_rnorm <- mock(c(0.6578201, 1.1364790), cycle = TRUE)
  with_mock(rnorm = mock_rnorm, {
    expect_equal(gaussian_system$get_momentum_sample(), L %*% c(0.6578201, 1.1364790))
  })
})
