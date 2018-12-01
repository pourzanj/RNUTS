context("Build Tree")

########################
# set up funnel system and inital point
########################
M <- diag(c(80.0, 9.0))
funnel_system <- create_funnel_hamiltonian_system(M, 1)

q0 <- c(-5.16181924428648, 0.13642758287681)
p0 <- c(-2.99685485641713, -0.14442408678326)
h0 <- 0.2

z0 <- list(q = q0, p = p0, h = h0)
H0 <- funnel_system$compute_H(z0)

########################
# test creation of a single node tree in isolation
########################
test_that("create_onenode_tree as a root tree has all the correct elements and it's history tibble is correct", {
  H <- funnel_system$compute_H(z0)

  onenode_tree <- create_onenode_tree(z = z0,
                              depth = NA,
                              H0 = H,
                              H = H,
                              valid_subtree = as.logical(NA),
                              uturn = as.logical(NA),
                              integrator_error = as.character(NA),
                              num_grad = as.integer(NA),
                              num_hess = as.integer(NA),
                              num_hess_vec = as.integer(NA),
                              num_newton = as.integer(NA),
                              DEBUG = TRUE)

  expect_true(is.na(onenode_tree$depth))
  expect_true(is.na(onenode_tree$valid))
  expect_true(is.null(onenode_tree$z_minus_1))
  expect_true(is.null(onenode_tree$z_minus_2))
  expect_true(is.null(onenode_tree$z_plus_1))
  expect_true(is.null(onenode_tree$z_plus_2))
})

test_that("create_onenode_tree works for the results of an explicit integrator", {
  take_leapfrog_step <- create_integrator(is_implicit = FALSE)
  integrator_result <- take_leapfrog_step(z0, z_1 = NULL, z_2 = NULL, direction = 1, funnel_system, H0)
  onenode_tree <- create_onenode_tree(z = integrator_result$z1,
                                  depth = 0,
                                  H0 = H0,
                                  H = funnel_system$compute_H(integrator_result$z1),
                                  valid_subtree = is.na(integrator_result$integrator_error),
                                  uturn = as.logical(NA),
                                  integrator_error = integrator_result$integrator_error,
                                  num_grad = integrator_result$num_grad,
                                  num_hess = integrator_result$num_hess,
                                  num_hess_vec = integrator_result$num_hess_vec,
                                  num_newton = integrator_result$num_newton,
                                  DEBUG = TRUE)

  expect_equal(onenode_tree$depth, 0)
  expect_true(onenode_tree$valid)
})

#
# test_that("create_onenode_tree works correctly when called by an implicit integrator", {
#   onenode_tree <- create_onenode_tree(depth = NA, invalid = NA, z = z0,
#                                       ham_system = funnel_system,
#                                       num_newton = NA, num_grad = NA, DEBUG = TRUE)
#
#   print(onenode_tree)
# })

########################
# test build_tree()
########################
test_that("base case of build_tree returns a single node tree", {

})
