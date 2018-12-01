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
create_integrator <- function(is_implicit, do_nonlin_solve = NULL) {

  # set integrator to either explicit or implicit
  # if implicit then the function we return
  integrate_one_step <- take_one_step_lf
  if (is_implicit) {
    integrate_one_step <- function(z0, z_1, z_2, direction, ham_system) {
      take_one_step_im(z0, z_1, z_2, direction, ham_system, do_nonlin_solve)
    }
  }

  # define and return stepper function
  function(z0, z_1, z_2, direction, ham_system, H0) {

    result <- integrate_one_step(z0, z_1, z_2, direction, ham_system, H0)

    list(z1 = result$z1,
         integrator_error = result$error,
         num_grad = result$num_grad,
         num_hess = ifelse(is_implicit, result$num_hess, as.integer(NA)),
         num_hess_vec = ifelse(is_implicit, result$num_hess, as.integer(NA)),
         num_newton = ifelse(is_implicit, result$num_newton, as.integer(NA)))
  }

}
