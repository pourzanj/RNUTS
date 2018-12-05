#' Single step implicit midpoint.
#'
#' Given a starting point, takes a single step using implicit midpoint.
#'
#' This function requires a function argument called do_nonlin_solve to solve
#' the nonlinear system that is required of implicit midpoint. do_nonlin_solve
#' can be as simple as a vanilla Newton iteration or more elaborate
#'
#' @param z0 starting point
#' @param z_1 previous point
#' @param z_2 previous, previous point
#' @param direction which direction to integrate in
#' @param ham_system the system we are integrating
#'
#' @return solution along with any errors and number of evaluations
#' @export
#'
#' @examples
take_one_step_im <- function(z0, z_1 = NULL, z_2 = NULL, direction, ham_system, do_nonlin_solve) {

  # unpack z0
  q0 <- z0$q
  p0 <- z0$p
  h <- z0$h*direction
  D <- length(q0) # dimension of q

  # unpack functions from ham_system
  M <- ham_system$M
  GradU <- function(q) ham_system$compute_gradU(q)
  HessU <- function(q) ham_system$compute_hessU(q)
  HessU_vec_prod <- function(q, v) ham_system$compute_hessU_vec_prod(q, v)

  # setup nonlinear solve functions required for backward euler half step
  # this includes function, its jacobian funciton, and a jacobian-vector product function
  g <- function(q_half) q_half - q0 - (h/2)*solve(M, p0 - (h/2)*GradU(q_half))
  Jg <- function(q_half) diag(D) + (h^2/4)*solve(M, HessU(q_half))
  Jg_v <- function(q_half, v) v + (h^2/4)*solve(M,HessU_vec_prod(q_half,v))

  # if we have two prev soln points use it to form a better Newton init guess
  x0 <- q0
  if (!is.null(z_1)) {
    q_1 <- z_1$q
    x0 <- ifelse(abs(q0+q_1)/abs(q0) < 0.1, 0.0, q0)
  }

  # do solve and extract info about num grad evals taken, whether newtons converged, etc
  nonlin_solve_soln <- do_nonlin_solve(g, Jg, Jg_v, x0)

  num_grad <- nonlin_solve_soln$num_func_evals
  num_hess <- nonlin_solve_soln$num_hess_evals
  num_hess_vec_prod <- nonlin_solve_soln$num_hess_vec_prod_evals
  num_newton <- nonlin_solve_soln$num_newton_iters
  error <- ifelse(nonlin_solve_soln$did_converge, NA, "Newton")

  # use solve to set backward euler half step
  q_half <- nonlin_solve_soln$x
  p_half <- p0 - (h/2)*GradU(q_half)

  # do forward euler half step
  q1 <- q.half + (h/2)*solve(M, p_half)
  p1 <- p.half - (h/2)*GradU(q_half)

  # pack back in to z1 and return
  z1 <- list(q = q1, p = p1, h = z0$h)

  list(z1 = z1, error = error, num_grad = num_grad, num_newton = num_newton,
       num_hess = num_hess, num_hess_vec_prod = num_hess_vec_prod)
}

#' Single step implicit midpoint in the midpoint formulation.
#'
#' Given a starting point, takes a single step using implicit midpoint. The solve is in the midpoint formulation rather than the backward euler, forward euler formulation.
#'
#' This function requires a function argument called do_nonlin_solve to solve
#' the nonlinear system that is required of implicit midpoint. do_nonlin_solve
#' can be as simple as a vanilla Newton iteration or more elaborate
#'
#' @param z0 starting point
#' @param z_1 previous point
#' @param z_2 previous, previous point
#' @param direction which direction to integrate in
#' @param ham_system the system we are integrating
#'
#' @return solution along with any errors and number of evaluations
#' @export
#'
#' @examples
take_one_step_im_midpoint <- function(z0, z_1 = NULL, z_2 = NULL, direction, ham_system, do_nonlin_solve) {

  # unpack z0
  q0 <- z0$q
  p0 <- z0$p
  h <- z0$h*direction
  D <- length(q0) # dimension of q

  # unpack functions from ham_system
  M <- ham_system$M
  GradU <- function(q) ham_system$compute_gradU(q)
  HessU <- function(q) ham_system$compute_hessU(q)
  HessU_vec_prod <- function(q, v) ham_system$compute_hessU_vec_prod(q, v)

  # setup nonlinear solve functions required for backward euler half step
  # this includes function, its jacobian funciton, and a jacobian-vector product function
  g <- function(p1) p1 - p0 + h*GradU(0.5*(q0 + q0 + h*solve(M, 0.5*(p0 + p1))))
  Jg <- function(p1) diag(D) + h*HessU(0.5*(q0 + q0 + h*solve(M, 0.5*(p0 + p1)))) %*% (h/4)*solve(M)
  Jg_v <- function(p1, v) v + h*HessU_vec_prod(0.5*(q0 + q0 + h*solve(M, 0.5*(p0 + p1))), (h/4)*solve(M, v))

  # use previous, previous point as initial guess. When the stepsize is bigger than the natural frequency
  # of the system, the numerical solution will oscillate so not the previous point, but the previous previous
  # point would be close. The previous previous point also however works when the solution is slow moving
  # because in that case the previous previous point is still a good guess
  x0 <- p0
  if (!is.null(z_2)) {
    p_2 <- z_2$p
    x0 <- p_2
  }

  # do solve and extract info about num grad evals taken, whether newtons converged, etc
  nonlin_solve_soln <- do_nonlin_solve(g, Jg, Jg_v, x0)

  num_grad <- nonlin_solve_soln$num_func_evals
  num_hess <- nonlin_solve_soln$num_hess_evals
  num_hess_vec_prod <- nonlin_solve_soln$num_hess_vec_prod_evals
  num_newton <- nonlin_solve_soln$num_newton_iters
  error <- ifelse(nonlin_solve_soln$did_converge, NA, "Newton")
  if(!nonlin_solve_soln$did_converge) {
    print("Newton did not converge")
  }

  # use solve to set q
  p1 <- nonlin_solve_soln$x
  q1 <- q0 + h*solve(M, 0.5*(p0 + p1))

  # pack back in to z1 and return
  z1 <- list(q = q1, p = p1, h = z0$h)

  list(z1 = z1, error = error, num_grad = num_grad, num_newton = num_newton,
       num_hess = num_hess, num_hess_vec_prod = num_hess_vec_prod)
}
