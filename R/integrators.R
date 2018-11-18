#' Title
#'
#' @param is_implicit
#' @param do_nonlin_solve a function that solves the nonlin system g(x) = 0
#'
#' @return
#' @export
#'
#' @examples
create_integrator <- function(is_implicit, do_nonlin_solve = NULL) {

  # set integrator to either explicit or implicit
  integrate_one_step <- take_one_step_lf
  if (is_implicit) {
    integrate_one_step <- function(z0, z_1, direction, ham_system) {
      take_one_step_im(z0, z_1, direction, ham_system, do_nonlin_solve)
    }
  }

  # define and return step function
  take_step <- function(z0, z_1, direction, ham_system, DEBUG) {

    result <- integrate_one_step(z0, z_1, direction, ham_system)
    H <- ham_system$compute_H(result$z1)

    num_newton <- as.integer(NA)
    num_grad <- result$num_grad
    if (is_implicit) {
      num_newton <- result$num_newton
    }

    create_onenode_tree(depth = 0,
                        invalid = result$error,
                        z = result$z1,
                        ham_system = ham_system,
                        num_newton = num_newton,
                        num_grad = num_grad,
                        DEBUG = DEBUG)
  }

  list(take_step = take_step)
}
