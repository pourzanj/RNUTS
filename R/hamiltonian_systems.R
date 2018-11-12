create_guassian_hamiltonian_system <- function(M, Sigma) {

  list(M = M,
       compute_H = function(z) (1/2)*sum(z$p*(1/M)*z$p) + compute_U(z$q),
       compute_gradU = compute_gradU,
       compute_hessU_vec_prod = compute_hessU_vec_prod,
       get_momentum_sample = function() L %*% rnorm(D))

}

create_funnel_hamiltonian_system <- function(M, D) {

}

#' Title
#'
#' @param M
#' @param compute_U
#' @param compute_gradU
#' @param compute_hessU_vec_prod
#'
#' @return
#' @export
#'
#' @examples
create_custom_hamiltonian_system <- function(M, compute_U, compute_gradU, compute_hessU_vec_prod) {

  if (is.matrix(M)) {
    stop("Only diaganol mass matrices specified as vectors are supported as of now.")
  }

  D <- length(M)
  L <- chol(diag(M))

  list(M = M,
       compute_H = function(z) (1/2)*sum(z$p*(1/M)*z$p) + compute_U(z$q),
       compute_gradU = compute_gradU,
       compute_hessU_vec_prod = compute_hessU_vec_prod,
       get_momentum_sample = function() L %*% rnorm(D))

}
