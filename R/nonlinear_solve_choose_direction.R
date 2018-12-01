#' Title
#'
#' @param g_x
#' @param Jg
#' @param Jg_v
#' @param x the point we're evaluating at
#'
#' @return
#' @export
#'
#' @examples
select_classic_newton_direction <- function(iter, x, g_x, Jg, Jg_v, compute_residual) {
  delta <- solve(Jg(x), -g_x)
  return(list(delta = delta, num_hess_evals = 1, num_hess_vec_prod_evals = 0))
}

#' Create a GMRES based direction_select function.
#'
#' A function to create a select_direction function based on GMRES that
#' uses a user-specified select_eta function e.g. one of the Eisenstat
#' Walker eta (linear solve tolerance) selectors or something more primitive.
#'
#' The rationale for gmres_x0_use_prev, is that delta is essentially a direction that we are updating in
#' so it might make sense that for several moves in a row delta can be in the same direction, in which
#' case if the initial delta guess delta0 is close to the new delta solution then the GMRES will have to
#' do fewer iterations in theory to get to the same low residual tolerance.
#'
#' Note that this function factory will create a select_direction() function once that will be used
#' for solving more than one nonlinear system, thus we must pass it the iteration number we're on
#' to keep track of whether we're on the 0th iteration and should just use the default eta_0.
#'
#' @param select_eta
#' @param gmres_x0_use_prev whether to use the previous nonlin soln as guess to GMRES or zero vector
#'
#' @return
#' @export
#'
#' @examples
create_gmres_direction_selector <- function(select_eta, gmres_delta0_use_prev = FALSE) {

  # closure variables thatwill get updated at each function call of select_direction
  # this is so that the function "remembers" what the previous values are so that it
  # can use it to select better values of eta ala Eisenstat and Walker

  # previous tolerance used
  eta_prev <- NULL

  # previous solution
  x_prev <- NULL

  # the nonlinear system evaluated at the previous point
  # starts as NULL but gets set in every function call
  g_x_prev <- NULL

  # delta_{k-1}
  delta_prev <- NULL

  # norm of the previous linear solution
  # (see equations 2.1 and 2.2 of Eisenstat and Walker)
  norm_lin_soln_res_previous <- NULL

  function(iter, x, g_x, Jg, Jg_v, compute_residual) {

    # use zero vector as initial guess to GMRES if gmres_delta0_use_prev is on
    # otherwise use previous solution
    delta0 <- rep(0.0, length(g_x))
    if(gmres_delta0_use_prev & iter > 0) {
      delta0 <- delta_prev
    }

    # keep track of how many function evaluations the select eta and actual GMRES solve are doing
    num_hess_vec_prod_evals <- 0

    # select eta which probably depends on closure variables
    selected_eta <- select_eta(iter, eta_prev, x, x_prev, g_x, g_x_prev, Jg_v, delta_prev, norm_lin_soln_res_previous)
    eta_prev <<- selected_eta$eta
    num_hess_vec_prod_evals <- num_hess_vec_prod_evals + selected_eta$num_hess_vec_prod_evals

    # get GMRES soln and extra info such as number of evals and tolerance
    # note that the gmres() function requires a matrix product function so we
    # make a small helper function that evaluations Jg at x then does the multiply
    Jg_x_v <- function(v) Jg_v(x, v)
    gmres_soln <- gmres(Jg_x_v, -g_x, delta0, TOL = selected_eta$eta)
    num_hess_vec_prod_evals <- num_hess_vec_prod_evals + gmres_soln$iterations

    # Update closure variables which are necessary for select_eta() on the next iteration of the Newton.
    # Note that Eigen GMRES's error is already
    # divided by the norm of b in Ax=b so we must multiply to get |Ax-b|
    # which is what Eisenstat 2.1 uses
    x_prev <<- x
    g_x_prev <<- g_x
    delta_prev <<- gmres_soln$x
    norm_lin_soln_res_previous <<- gmres_soln$error * compute_two_norm(g_x)

    return(list(delta = gmres_soln$x, num_func_evals = 0, num_hess_evals = 0, num_hess_vec_prod_evals = num_hess_vec_prod_evals))
  }

}

select_eta_constant_1e1 <- function(iter, eta_prev, x, x_prev, g_x, g_x_prev, Jg_v, delta_prev, norm_lin_soln_res_previous) {
  list(eta = 1e-1, num_hess_vec_prod_evals = 0)
}

select_eta_constant_1e4 <- function(iter, eta_prev, x, x_prev, g_x, g_x_prev, Jg_v, delta_prev, norm_lin_soln_res_previous) {
  list(eta = 1e-4, num_hess_vec_prod_evals = 0)
}

select_eta_saad <- function(iter, eta_prev, x, x_prev, g_x, g_x_prev, Jg_v, delta_prev, norm_lin_soln_res_previous) {
  list(eta = 1/2^(iter+1), num_hess_vec_prod_evals = 0)
}

#' Create Eisenstat and Walker Choice 1
#'
#' Creates the choice 1 select_eta function from Eisenstat and Walker. Using
#' either choice 2.1 or 2.2.
#'
#' @param choice Whether to use 2.1 or 2.2 from Eisenstat and Walker
#'
#' @return
#' @export
#'
#' @examples
create_select_eta_choice_1 <- function(choice = "2.2") {

  function(iter, eta_prev, x, x_prev, g_x, g_x_prev, Jg_v, delta_prev, norm_lin_soln_res_previous) {

    # keep track of function evaluations (actually just for choice 2.1)
    num_hess_vec_prod_evals <- 0

    # eta will always start at 0.5 for the 0th iteration. per Eisenstat and Walker
    eta <- 0.5

    # for subsequent iterations run the selector either choice 2.1 or 2.2
    if(iter > 0) {
      if(choice == "2.1") {
        eta <- compute_two_norm(g_x - g_x_prev - Jg_v(x_prev, delta_prev)) / compute_two_norm(g_x_prev)
        num_hess_vec_prod_evals <- num_hess_vec_prod_evals + 1
      }
      else {
        eta <- abs(compute_two_norm(g_x)-norm_lin_soln_res_previous) / compute_two_norm(g_x_prev)
      }

      # use safeguard to make sure forcing term doesn't get small to quickly
      # but only near the beginning of the nonlinear solve where the etas are still above 0.1
      r <- (1.0+sqrt(5.0))/2.0
      if(eta_prev^r > 0.1) {
        eta <- max(eta, eta_prev^r)
      }
    }

    #print(paste0("eta: ", eta))

    list(eta = eta, num_hess_vec_prod_evals = num_hess_vec_prod_evals)
  }

}

#' Create select_eta choice 2 from Eisentstat and Walker with a given set of hyperparameters.
#'
#' Unlike the other select_eta functions, choice 2 has hyperparameters that we have to set in the beginning.
#'
#' @return
#' @export
#'
#' @examples
select_eta_choice_2 <- function(gamma = 0.9, alpha = (1.0+sqrt(5.0))/2.0) {

  function(iter, eta_prev, x, x_prev, g_x, g_x_prev, Jg_v, delta_prev, norm_lin_soln_res_previouss) {

    # eta will always start at 0.5 for the 0th iteration. per Eisenstat and Walker
    eta <- 0.5
    if(iter > 0) {
      eta <- gamma*(compute_two_norm(g_x)/compute_two_norm(g_x_prev))^alpha

      # use safeguard to make sure forcing term doesn't get small to quickly
      # but only near the beginning of the nonlinear solve where the etas are still above 0.1
      if(gamma*eta_prev^alpha > 0.1) {
        eta <- max(eta, gamma*eta_prev^alpha)
      }
    }

    list(eta = eta, num_hess_vec_prod_evals = 0)
  }
}

