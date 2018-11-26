create_constant_select_stepsize <- function(alpha = 1.0) {
  function(iter, delta, g_x, g, Jg, Jg_v, x, compute_residual) {
    list(alpha = alpha, num_func_evals = 0, num_hess_evals = 0, num_hess_vec_prod_evals = 0)
  }
}

create_geometric_select_stepsize <- function(gamma = 2.0) {
  function(iter, delta, g_x, g, Jg, Jg_v, x, compute_residual) {
    list(alpha = gamma^iter, num_func_evals = 0, num_hess_evals = 0, num_hess_vec_prod_evals = 0)
  }
}


