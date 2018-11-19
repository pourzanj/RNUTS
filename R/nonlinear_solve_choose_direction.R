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
select_classic_newton_direction <- function(iter, g_x, Jg, Jg_v, x, compute_residual) {
  delta <- solve(Jg(x), -g_x)
  return(list(delta = delta, num_hess_evals = 1, num_hess_vec_prod_evals = 0))
}

#' Create a GMRES based direction_select function.
#'
#' A function to create a select_direction function based on GMRES that
#' uses a user-specified select_eta function e.g. one of the Eisenstat
#' Walker eta (linear solve tolerance) selectors or something more primitive.
#'
#' The rational for gmres_x0_use_prev, is that delta is essentially a direction that we are updating in
#' so it might make sense that for several moves in a row delta can be in the same direction, in which
#' case if the initial delta guess delta0 is close to the new delta solution then the GMRES will have to
#' do fewer iterations in theory to get to the same low residual tolerance.
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

  # the nonlinear system evaluated at the previous point
  # starts as NULL but gets set in every function call
  g_x_prev <- NULL

  # previous linear solution the direction select chose
  delta_prev <- NULL

  # Jg(x_{k-1}) %*% delta_{k-1}
  Jg_prev_delta_prev <- NULL

  # norm of the previous linear solution
  # (see equations 2.1 and 2.2 of Eisenstat and Walker)
  norm_lin_soln_res_previous <- NULL

  function(iter, g_x, Jg, Jg_v, x, compute_residual) {

    # use zero vector based if gmres_delta0_use_prev is on
    delta0 <- rep(0.0, length(g_x))
    if(gmres_delta0_use_prev) {
      delta0 <- delta_prev
    }

    # select eta based on closure variables and given select_etc function
    # then use output to update closure variables
    eta <- select_eta(iter, g_x, g_x_prev, delta_prev, Jg_prev_delta_prev, norm_lin_soln_res_previous, compute_residual)
    g_x_prev <<- eta$g_x_prev
    delta_prev <<- eta$delta_prev
    Jg_prev_delta_prev <<- eta$Jg_prev_delta_prev

    # get GMRES soln and extra info such as number of evals and tolerance
    # the error could be useful for Eisentstat and Walker eq. 2.2
    gmres_soln <- gmres(Jg_v, -g_x, x0, TOL = eta$eta)
    norm_lin_soln_res_previous <<- gmres_soln$error

    return(list(delta = gmres_soln$x, num_hess_evals = 0, num_hess_vec_prod_evals = gmres_soln$iterations))
  }

}

select_eta_constant <- function(iter, g_x, g_x_prev, delta_prev, Jg_prev_delta_prev, norm_lin_soln_res_previous, compute_residual) {
  list(eta = eta, g_x_prev = NULL, delta_prev = NULL, Jg_prev_delta_prev = NULL)
}

select_eta_saad <- function(iter, g_x, g_x_prev, delta_prev, Jg_prev_delta_prev, norm_lin_soln_res_previous, compute_residual) {
  list(eta = eta, g_x_prev = NULL, delta_prev = NULL, Jg_prev_delta_prev = NULL)
}

select_eta_choice_1A <- function(iter, g_x, g_x_prev, delta_prev, Jg_prev_delta_prev, norm_lin_soln_res_previous, compute_residual) {
  list(eta = eta, g_x_prev = NULL, delta_prev = NULL, Jg_prev_delta_prev = NULL)
}

select_eta_choice_1B <- function(iter, g_x, g_x_prev, delta_prev, Jg_prev_delta_prev, norm_lin_soln_res_previous, compute_residual) {
  list(eta = eta, g_x_prev = NULL, delta_prev = NULL, Jg_prev_delta_prev = NULL)
}

#' Create select_eta choice 2 from Eisentstat and Walker with a given set of hyperparameters.
#'
#' Unlike the other select_eta functions, choice 2 has hyperparameters that we have to set in the beginning.
#'
#' @return
#' @export
#'
#' @examples
select_eta_choice_2 <- function() {
  function(iter, g_x, g_x_prev, delta_prev, Jg_prev_delta_prev, norm_lin_soln_res_previous, compute_residual) {
    list(eta = eta, g_x_prev = NULL, delta_prev = NULL, Jg_prev_delta_prev = NULL)
  }
}

