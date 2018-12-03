context("Creating an Integrator")

M <- diag(c(9.0, 80.0))
funnel_system <- create_funnel_hamiltonian_system(M, 1)

q0 <- c(-5.16181924428648, 0.13642758287681)
p0 <- c(-2.99685485641713, -0.14442408678326)
h0 <- 0.5

take_lf_step <- create_integrator(is_implicit = FALSE)

test_that("integrate_multiple_steps() has results that are correctly formatted", {
  soln <- integrate_multiple_steps(q0, p0, h0, funnel_system, take_lf_step, num_steps = 10)
  expect_equal(soln$num_grad[2], 2)
})

test_that("integrate_multiple_steps() stops and reports divergences ", {
  h0 <- 5.0
  soln <- integrate_multiple_steps(q0, p0, h0, funnel_system, take_lf_step, num_steps = 10)
  expect_equal(soln$integrator_error[2], "Divergence")
})
