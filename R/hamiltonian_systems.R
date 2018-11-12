create_gaussian_hamiltonian_system <- function(M, Sigma) {

  compute_U <- function(x) 0.5*sum(x*solve(Sigma,x))
  compute_gradU <- function(x) solve(Sigma,x)
  compute_hessU_vec_prod <- function(q, v) solve(Sigma, v)
  D <- nrow(M)
  L <- chol(M)

  list(M = M,
       compute_H = function(z) 0.5*sum(z$p*solve(M,z$p)) + compute_U(z$q),
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
  L <- chol(M)

  list(M = M,
       compute_H = function(z) (1/2)*sum(z$p*(1/M)*z$p) + compute_U(z$q),
       compute_gradU = compute_gradU,
       compute_hessU_vec_prod = compute_hessU_vec_prod,
       get_momentum_sample = function() L %*% rnorm(D))

}
