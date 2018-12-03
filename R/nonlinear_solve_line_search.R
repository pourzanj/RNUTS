create_constant_select_stepsize <- function(alpha = 1.0) {
  function(iter, delta, g_x, g, Jg, Jg_v, x) {
    list(alpha = alpha, num_func_evals = 0, num_hess_evals = 0, num_hess_vec_prod_evals = 0)
  }
}

create_geometric_select_stepsize <- function(gamma = 0.5, MAX_ITER = 10) {
  function(iter, delta, g_x, g, Jg, Jg_v, x) {

    alpha <- 1.0
    num_func_evals <- 0

    g_x_new <- g(x + alpha*delta)
    norm_g_x <- compute_two_norm(g_x)
    norm_g_x_new <- compute_two_norm(g_x_new)
    num_func_evals <- num_func_evals + 1

    while(norm_g_x_new > norm_g_x) {
      alpha <- gamma*alpha

      g_x_new <- g(x + alpha*delta)
      norm_g_x_new <- compute_two_norm(g_x_new)
      num_func_evals <- num_func_evals + 1
    }


    list(alpha = alpha, num_func_evals = num_func_evals, num_hess_evals = 0, num_hess_vec_prod_evals = 0)
  }
}
