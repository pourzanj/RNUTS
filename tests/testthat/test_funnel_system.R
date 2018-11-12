library(mockery)

context("The Hamiltonian System for a Funnel")

######################
# Simple 1D Funnel
######################
M <- diag(c(80.0, 9.0))
funnel_system <- create_funnel_hamiltonian_system(M, 1)

z1 <- list(q = c(0.0, 0.0), p = c(-0.1759228, 0.8358274))
z2 <- list(q = c(0.08784589, -0.96785614), p = c(-0.1759228, 0.8358274))

test_that("M is included and matches input", {
  expect_equal(funnel_system$M, M)
})

test_that("Hamiltonian is computed correctly", {
  expect_equal(funnel_system$compute_H(z1), 0.5*sum(z1$p*solve(M,z1$p)))
  expect_equal(funnel_system$compute_H(z2), 0.5*sum(z2$p*solve(M,z2$p)) + 0.4733352)
})

test_that("gradU is computed correctly", {
  expect_equal(funnel_system$compute_gradU(z1$q), c(0.5, 0.0))
  expect_equal(funnel_system$compute_gradU(z2$q), c(0.08077711, -0.88646139))
})

test_that("hessU vector-product is computed correctly", {
  expect_equal(funnel_system$compute_hessU_vec_prod(z1$q, c(0.0, 0.0)), c(0.0, 0.0))
  expect_equal(funnel_system$compute_hessU_vec_prod(z2$q, c(0.5, 0.123)), c(0.37908208108335, 0.55588664258090))
})

test_that("momentum sampler is included and draws from the correct momentum", {
  L <- chol(M)
  mock_rnorm <- mock(c(0.6578201, 1.1364790), cycle = TRUE)
  with_mock(rnorm = mock_rnorm, {
    expect_equal(funnel_system$get_momentum_sample(), L %*% c(0.6578201, 1.1364790))
  })
})
