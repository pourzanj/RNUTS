library(mockery)

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
                              valid_subtree = TRUE,
                              uturn = FALSE,
                              integrator_error = as.character(NA),
                              num_grad = as.integer(NA),
                              num_hess = as.integer(NA),
                              num_hess_vec = as.integer(NA),
                              num_newton = as.integer(NA),
                              DEBUG = TRUE)

  expect_true(is.na(onenode_tree$depth))
  expect_true(onenode_tree$valid)
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
                                  uturn = FALSE,
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
# test build_tree() base case
########################
test_that("base case of build_tree returns a single node tree", {

  take_leapfrog_step <- create_integrator(is_implicit = FALSE)
  new_onenode_tree <- build_tree(depth = 0, z0, z_1 = NULL, z_2 = NULL, direction = 1, funnel_system, H0, integrate_step = take_leapfrog_step, DEBUG = TRUE)

  expect_equal(new_onenode_tree$depth, 0)

  # for a one node tree z_minus and z_plus and z_rep are all the same
  expect_identical(new_onenode_tree$z_minus, new_onenode_tree$z_plus)
})

########################
# test join_subtrees() on simple case of joining two onenode_trees
########################
test_that("we can properly sample a representative from two newly created onenode_trees", {
    # create root node
    H <- funnel_system$compute_H(z0)

    root_tree <- create_onenode_tree(z = z0,
                                        depth = NA,
                                        H0 = H,
                                        H = H,
                                        valid_subtree = TRUE,
                                        uturn = FALSE,
                                        integrator_error = as.character(NA),
                                        num_grad = as.integer(NA),
                                        num_hess = as.integer(NA),
                                        num_hess_vec = as.integer(NA),
                                        num_newton = as.integer(NA),
                                        DEBUG = TRUE)


    # create new onenode_tree from a step
    take_leapfrog_step <- create_integrator(is_implicit = FALSE)
    new_onenode_tree <- build_tree(depth = 0, z0, z_1 = NULL, z_2 = NULL, direction = 1, funnel_system, H0, integrate_step = take_leapfrog_step, DEBUG = TRUE)

    # the new step has less energy so we should sample it w/ prob. 1
    new_z_rep <- sample_new_representative(root_tree, new_onenode_tree)

    expect_identical(new_onenode_tree$z_rep, new_z_rep)
})

test_that("we can properly join a root tree with a newly create onenode_tree", {

  # create root node
  H <- funnel_system$compute_H(z0)

  root_tree <- create_onenode_tree(z = z0,
                                      depth = NA,
                                      H0 = H,
                                      H = H,
                                      valid_subtree = TRUE,
                                      uturn = FALSE,
                                      integrator_error = as.character(NA),
                                      num_grad = as.integer(NA),
                                      num_hess = as.integer(NA),
                                      num_hess_vec = as.integer(NA),
                                      num_newton = as.integer(NA),
                                      DEBUG = TRUE)


  # create new onenode_tree from a step
  take_leapfrog_step <- create_integrator(is_implicit = FALSE)
  new_onenode_tree <- build_tree(depth = 0, z0, z_1 = NULL, z_2 = NULL, direction = 1, funnel_system, H0, integrate_step = take_leapfrog_step, DEBUG = TRUE)

  # join the two
  joined_tree <- join_subtrees(root_tree, new_onenode_tree, direction = 1, funnel_system, DEBUG = TRUE)

  # z_plus actually lost energy so it should be the represntative w. prob. 1
  expect_identical(joined_tree$z_rep, joined_tree$z_plus)

  expect_true(joined_tree$valid)
})

########################
# test build_tree() on simple case of joining two onenode_trees
########################
