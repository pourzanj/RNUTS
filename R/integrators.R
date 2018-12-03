#' Create an integrator.
#'
#' Create either a Leapfrog integrator or an Implicit Midpoint integrator with a particular nonlinear solve function.
#'
#' @param is_implicit Whether to use Implicit Midpoint or Leapfrog
#' @param do_nonlin_solve a function that solves the nonlin system g(x) = 0
#'
#' @return A function that integrates a single step.
#' @export
#'
#' @examples
create_integrator <- function(is_implicit, midpoint_formulation = TRUE, do_nonlin_solve = NULL) {

  # set integrator to either explicit or implicit
  # if implicit then the function we return
  integrate_one_step <- take_one_step_lf
  if (is_implicit) {
    # use either Backward Euler formulation or Midpoint Formulation
    if(midpoint_formulation) {
      integrate_one_step <- function(z0, z_1, z_2, direction, ham_system, H0) {
        take_one_step_im_midpoint(z0, z_1, z_2, direction, ham_system, do_nonlin_solve)
      }
    } else {
      integrate_one_step <- function(z0, z_1, z_2, direction, ham_system, H0) {
        take_one_step_im(z0, z_1, z_2, direction, ham_system, do_nonlin_solve)
      }
    }
  }

  # define and return stepper function
  function(z0, z_1, z_2, direction, ham_system, H0) {

    result <- integrate_one_step(z0, z_1, z_2, direction, ham_system, H0)

    list(z1 = result$z1,
         integrator_error = result$error,
         num_grad = result$num_grad,
         num_hess = ifelse(is_implicit, result$num_hess, as.integer(NA)),
         num_hess_vec = ifelse(is_implicit, result$num_hess_vec_prod, as.integer(NA)),
         num_newton = ifelse(is_implicit, result$num_newton, as.integer(NA)))
  }

}

#' Integrate multiple steps.
#'
#' Uses a take step function to integrate multiple steps then format the result as a tibble.
#'
#' @param q0 initial position
#' @param p0 initial momentum
#' @param h0 initial stepsize
#' @param ham_system the Hamiltonian system to integrate
#' @param take_step the take single step function
#' @param num_steps the number of steps to integrate
#'
#' @return a formatted tibble containing information of the steps
#' @export
#'
#' @examples
integrate_multiple_steps <- function(q0, p0, h0, ham_system, take_step, num_steps) {

  D <- length(q0)
  h <- rep(NA, num_steps+1)
  H <- rep(NA, num_steps+1)
  integrator_error <- rep(NA, num_steps+1)
  num_grad <- rep(NA, num_steps+1)
  num_hess <- rep(NA, num_steps+1)
  num_hess_vec <- rep(NA, num_steps+1)
  num_newton <- rep(NA, num_steps+1)
  q <- matrix(NA, nrow = num_steps+1, ncol = D)
  p <- matrix(NA, nrow = num_steps+1, ncol = D)

  h[1] <- h0
  H[1] <- ham_system$compute_H(list(q = q0, p = p0, h = h0))
  q[1,] <- q0
  p[1,] <- p0

  # start solving
  for(n in 1:num_steps) {

    # pack variables
    z0 <- list(q = q[n,], p = p[n,], h = h[n])
    z_1 <- NULL
    z_2 <- NULL
    if(n > 2) {
      z_1 <- list(q = q[n-1,], p = p[n-1,], h = h[n-1])
      z_2 <- list(q = q[n-2,], p = p[n-2,], h = h[n-2])
    }

    # take step and catch any integrator errors
    one_step_result <- take_step(z0, z_1, z_2, direction = 1, ham_system, H[1])

    # unpack result and save
    h[n+1] <- one_step_result$z1$h
    H[n+1] <- ham_system$compute_H(list(q = one_step_result$z1$q, p = one_step_result$z1$p, h = one_step_result$z1$h))
    integrator_error[n+1] <- one_step_result$integrator_error
    num_grad[n+1] <- one_step_result$num_grad
    num_hess[n+1] <- one_step_result$num_hess
    num_hess_vec[n+1] <- one_step_result$num_hess_vec
    num_newton[n+1] <- one_step_result$num_newton
    q[n+1,] <- one_step_result$z1$q
    p[n+1,] <- one_step_result$z1$p

    if(!is.na(one_step_result$integrator_error)) break
  }

  # return as a tibble
  soln <- tibble(t = cumsum(c(0,tail(h, num_steps)))) %>%
    mutate(h = h,
           H0 = H[1],
           H = H,
           integrator_error = integrator_error,
           num_grad = num_grad,
           num_hess = num_hess,
           num_hess_vec = num_hess_vec,
           num_newton = num_newton) %>%
    bind_cols(as_tibble(q) %>% set_names(paste0("q", 1:D))) %>%
    bind_cols(as_tibble(p) %>% set_names(paste0("p", 1:D)))

  return(soln)
}
