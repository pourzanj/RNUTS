context("Implicit Midpoint")

test_that("Implicit Midpoint in the Midpoint Formulation Steps Correctly on 2D Gaussian", {
  M <- diag(c(1.0, 1.0))
  Sigma <- matrix(c(1.0, 0.9, 0.9, 1.0), nrow = 2)
  gaussian_system <- create_gaussian_hamiltonian_system(M, Sigma)
  q0 <- c(0.0, 0.0)
  p0 <- c(-0.1759228, 0.8358274)
  h0 <- 0.5

  # set up nonlinear solve
  select_eta_choice_1_2_2 <- create_select_eta_choice_1(choice = "2.2")
  select_direction_gmres_choice_1_2_2 <- create_gmres_direction_selector(select_eta_choice_1_2_2 , gmres_delta0_use_prev = FALSE)
  select_stepsize_geometric_backtrack <- create_geometric_select_stepsize(gamma = 0.5, MAX_ITER = 10)

  newton_gmres_solve <- create_generalized_newton_solve(select_direction_gmres_choice_1_2_2,
                                                        select_stepsize_geometric_backtrack,
                                                        compute_rms_norm,
                                                        NTOL = 1e-6, MAX_ITER = 20)

  # set up integrator
  take_im_midpoint_step <- create_integrator(is_implicit = TRUE, midpoint_formulation = TRUE, do_nonlin_solve = newton_gmres_solve)

  # get soln for multiple steps
  soln <- integrate_multiple_steps(q0, p0, h0, gaussian_system, take_im_midpoint_step, num_steps = 10)

  expect_true(!is.na(soln$num_hess[2]))

  # plot
  soln %>%
    select(t, H, q1, q2, p1, p2) %>%
    gather(state, value, -t) %>%
    ggplot(aes(t, value)) +
    geom_point() +
    geom_line() +
    facet_grid(state ~ ., scales = "free")
})

test_that("Implicit Midpoint in the Midpoint Formulation Steps Correctly on 2D Funnel", {

  M <- diag(c(9.0, 80.0))
  funnel_system <- create_funnel_hamiltonian_system(M, 1)
  q0 <- c(-5.16181924428648, 0.13642758287681)
  p0 <- c(-2.99685485641713, -0.14442408678326)
  h0 <- 3.0

  # set up nonlinear solve
  select_eta_choice_1_2_2 <- create_select_eta_choice_1(choice = "2.2")
  select_direction_gmres_choice_1_2_2 <- create_gmres_direction_selector(select_eta_choice_1_2_2 , gmres_delta0_use_prev = FALSE)
  select_stepsize_geometric_backtrack <- create_geometric_select_stepsize(gamma = 0.5, MAX_ITER = 10)

  newton_gmres_solve <- create_generalized_newton_solve(select_direction_gmres_choice_1_2_2,
                                                        select_stepsize_geometric_backtrack,
                                                        compute_rms_norm,
                                                        NTOL = 1e-12, MAX_ITER = 20)

  # set up integrator
  take_im_midpoint_step <- create_integrator(is_implicit = TRUE, midpoint_formulation = TRUE, do_nonlin_solve = newton_gmres_solve)

  # get soln for multiple steps
  soln <- integrate_multiple_steps(q0, p0, h0, funnel_system, take_im_midpoint_step, num_steps = 10)

  expect_true(!is.na(soln$num_hess[2]))

  # plot
  soln %>%
    select(t, H, q1, q2, p1, p2) %>%
    gather(state, value, -t) %>%
    ggplot(aes(t, value)) +
    geom_point() +
    geom_line() +
    facet_grid(state ~ ., scales = "free")

  # check that it is reversible i.e. that we get back to the same point when we integrate backwards
  q0_backwards <- soln %>% select(q1, q2) %>% tail(1) %>% as.matrix() %>% as.vector()
  p0_backwards <- soln %>% select(p1, p2) %>% tail(1) %>% as.matrix() %>% as.vector()
  soln_backwards <- integrate_multiple_steps(q0_backwards, -p0_backwards, h0, funnel_system, take_im_midpoint_step, num_steps = 10)
  q_backwards <- soln_backwards %>% select(q1, q2) %>% tail(1) %>% as.matrix() %>% as.vector()
  p_backwards <- soln_backwards %>% select(p1, p2) %>% tail(1) %>% as.matrix() %>% as.vector()

  expect_equal(q_backwards, q0)
  expect_equal(-p_backwards, p0)
})

