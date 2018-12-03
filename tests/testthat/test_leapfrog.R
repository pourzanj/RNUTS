context("Leapfrog")

test_that("Onestep Leapfrog Steps Correctly on 2D Gaussian", {
  M <- diag(c(1.0, 1.0))
  Sigma <- matrix(c(1.0, 0.9, 0.9, 1.0), nrow = 2)
  gaussian_system <- create_gaussian_hamiltonian_system(M, Sigma)

  # compute Leapfrog step manually
  q0 <- c(0.0, 0.0)
  p0 <- c(-0.1759228, 0.8358274)
  p_half <- p0 - (0.1/2)*gaussian_system$compute_gradU(q0)
  q1 <- q0 + 0.1*solve(M, p_half)
  p1 <- p_half - (0.1/2)*gaussian_system$compute_gradU(q1)

  # compute using library function
  z0 <- list(q = c(0.0, 0.0), p = c(-0.1759228, 0.8358274), h = 0.1)
  H0 <- gaussian_system$compute_H(z0)
  result <- take_one_step_lf(z0, z_1 = NULL, z_2 = NULL, direction = 1, gaussian_system, H0)
  z1 <- result$z1

  expect_equal(z1$q, q1)
  expect_equal(z1$p, p1)
  expect_equal(z1$h, 0.1)
})

test_that("Onestep Leapfrog Catches Divergences", {
  M <- diag(c(1.0, 1.0))
  Sigma <- matrix(c(1.0, 0.9, 0.9, 1.0), nrow = 2)
  gaussian_system <- create_gaussian_hamiltonian_system(M, Sigma)

  z0 <- list(q = c(0.0, 0.0), p = c(-0.1759228, 0.8358274), h = 0.1)
  H0 <- gaussian_system$compute_H(z0)
  result <- take_one_step_lf(z0, z_1 = NULL, z_2 = NULL, direction = 1, gaussian_system, H0)
  z1 <- result$z1

  # try with a big step that should be divergent
  z0 <- list(q = c(0.0, 0.0), p = c(-0.1759228, 0.8358274), h = 10.0)
  H0 <- gaussian_system$compute_H(z0)

  result <- take_one_step_lf(z0, z_1 = NULL, z_2 = NULL, direction = 1, gaussian_system, H0)

  expect_true(!is.na(result$error))
  expect_identical(result$error, "Divergence")
})
