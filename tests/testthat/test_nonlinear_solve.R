context("Nonlinear Solve")

###############
# SETUP FUNNEL
###############

# set up a Funnel system to test a typical nonlinear solve that we'll
# have to do when using Implicit Midpoint
M <- diag(c(80.0, 9.0))
funnel_system <- create_funnel_hamiltonian_system(M, 1)
GradU <- function(q) funnel_system$compute_gradU(q)
HessU <- function(q) funnel_system$compute_hessU(q)
HessU_vec_prod <- function(q, v) funnel_system$compute_hessU_vec_prod(q, v)

q0 <- c(-5.16181924428648, 0.13642758287681)
p0 <- c(-2.99685485641713, -0.14442408678326)
h <- 0.1

D <- 2
g <- function(q_half) q_half - q0 - (h/2)*solve(M, p0 - (h/2)*GradU(q_half))
Jg <- function(q_half) diag(D) + (h^2/4)*solve(M, HessU(q_half))
Jg_v <- function(q_half, v) v + (h^2/4)*solve(M,HessU_vec_prod(q_half,v))


###############
# TESTS
###############
test_that("Generalized nonlinear solve works with most basic settings i.e. classic Newton and no line search", {

  select_stepsize_constant_one <- create_constant_select_stepsize(1.0)
  classic_newton_solve <- create_generalized_newton_solve(select_classic_newton_direction,
                                                   select_stepsize_constant_one,
                                                   compute_rms_norm,
                                                   NTOL = 1e-6, MAX_ITER = 20)

  nonlin_solve_soln <- classic_newton_solve(g, Jg, Jg_v, x0 = q0)
  eps <- compute_rms_norm(g(nonlin_solve_soln$x))
  expect_true(eps < 1e-6)
})
