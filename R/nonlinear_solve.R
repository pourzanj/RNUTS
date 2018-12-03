#' Create generalized Newton solve function.
#'
#' A function for creating a custom function that solves g(x)=0. Special cases include Newton iterations and Newton-Krylov.
#'
#' This function returns a function to solve nonlinear equations using a given select_direction,
#' and select_stepsize function. Internally, the function will be updating the intial guess as
#' x = x + alpha*delta.
#'
#' The select_direction function selects the vector delta. For Newton iterations
#' this would simply be delta = solve(Jg(x), -g(x)).
#'
#' The select_stepsize function selects alpha. For plain Newton iterations this would always be
#' one, but it can be replaced by more elaborate techniques such as line search.
#'
#' @param select_direction
#' @param select_stepsize
#' @param TOL
#' @param MAX_ITER
#'
#' @return
#' @export
#'
#' @examples
create_generalized_newton_solve <- function(select_direction, select_stepsize, compute_residual, NTOL = 1e-6, MAX_ITER = 20) {

  function(g, Jg, Jg_v, x0) {

    iter <- 0
    num_func_evals <- 0
    num_hess_evals <- 0
    num_hess_vec_prod_evals <- 0

    x <- x0
    g_x <- g(x)
    num_func_evals <- num_func_evals + 1
    eps <- compute_residual(g_x)
    converged <- eps < NTOL

    while(!converged) {

      #print(paste0("Iteration: ", iter, " eps: ", eps))

      # 1) choose direction to update, delta
      direction_soln <- select_direction(iter, x, g_x, Jg, Jg_v, compute_residual)
      delta <- direction_soln$delta
      num_hess_evals <- num_hess_evals + direction_soln$num_hess_evals
      num_hess_vec_prod_evals <- num_hess_vec_prod_evals + direction_soln$num_hess_vec_prod_evals

      # 2) choose stepsize
      stepsize_soln <- select_stepsize(iter, delta, g_x, g, Jg, Jg_v, x)
      alpha <- stepsize_soln$alpha
      num_func_evals <- num_func_evals + stepsize_soln$num_func_evals
      num_hess_evals <- num_hess_evals + stepsize_soln$num_hess_evals
      num_hess_vec_prod_evals <- num_hess_vec_prod_evals + stepsize_soln$num_hess_vec_prod_evals

      # 3) make update
      x <- x + alpha*delta

      # 4) compute new error and check if we've reached the max number of iterations
      g_x <- g(x)
      num_func_evals <- num_func_evals + 1

      eps <- compute_residual(alpha*delta)
      converged <- eps < NTOL

      iter <- iter + 1
      if(iter >= MAX_ITER & !converged) {
        warning("Newton reached max number of iterations and failed to converge.")
        break
      }
    }

    return(list(x = x, did_converge = converged, num_newton_iters = iter,
                num_func_evals = num_func_evals,
                num_hess_evals = num_hess_evals,
                num_hess_vec_prod_evals = num_hess_vec_prod_evals))
  }

}
