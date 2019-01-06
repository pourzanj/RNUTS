create_gaussian_hamiltonian_system <- function(M, Sigma) {

  compute_U <- function(x) 0.5*sum(x*solve(Sigma,x))
  compute_gradU <- function(x) solve(Sigma,x)
  compute_hessU <- function(x) solve(Sigma)
  compute_hessU_vec_prod <- function(q, v) solve(Sigma, v)
  D <- nrow(M)
  L <- chol(M)

  list(M = M,
       compute_U = compute_U,
       compute_H = function(z) 0.5*sum(z$p*solve(M,z$p)) + compute_U(z$q),
       compute_gradU = compute_gradU,
       compute_hessU = compute_hessU,
       compute_hessU_vec_prod = compute_hessU_vec_prod,
       get_momentum_sample = function() L %*% rnorm(D))

}

create_funnel_hamiltonian_system <- function(M, D) {

  compute_U <- function(q) {
    y <- q[1]
    x <- q[2:length(q)]
    N <- length(q) - 1
    (y^2)/18 + N*(y/2) + sum(x^2)/(2*exp(y))
  }

  compute_gradU <- function(q) {
    y <- q[1]
    x <- q[2:length(q)]
    N <- length(q) - 1
    c(y/9 + N/2 - sum(x^2)/(2*exp(y)), x*exp(-y))
  }

  compute_hessU <- function(q) {
    y <- q[1]
    x <- q[2:length(q)]

    # arrow matrix
    H <- diag(length(q))
    diag(H) <- c(1/9 + sum(x^2)/(2*exp(y)), rep(exp(-y), length(q)-1))
    H[2:length(q), 1] <- -x*exp(-y)
    H[1,2:length(q)] <- -x*exp(-y)

    H
  }

  compute_hessU_vec_prod <- function(q, v) {
    y <- q[1]
    x <- q[2:length(q)]

    # arrow matrix
    H <- diag(length(q))
    diag(H) <- c(1/9 + sum(x^2)/(2*exp(y)), rep(exp(-y), length(q)-1))
    H[2:length(q), 1] <- -x*exp(-y)
    H[1,2:length(q)] <- -x*exp(-y)

    as.vector(H %*% v)
  }

  D <- nrow(M)
  L <- chol(M)

  list(M = M,
       compute_U = compute_U,
       compute_H = function(z) 0.5*sum(z$p*solve(M,z$p)) + compute_U(z$q),
       compute_gradU = compute_gradU,
       compute_hessU = compute_hessU,
       compute_hessU_vec_prod = compute_hessU_vec_prod,
       get_momentum_sample = function() L %*% rnorm(D))
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
create_custom_hamiltonian_system <- function(M, compute_U, compute_gradU, compute_hessU, compute_hessU_vec_prod) {

  D <- nrow(M)
  L <- chol(M)

  list(M = M,
       compute_U = compute_U,
       compute_H = function(z) 0.5*sum(z$p*solve(M,z$p)) + compute_U(z$q),
       compute_gradU = compute_gradU,
       compute_hessU = compute_hessU,
       compute_hessU_vec_prod = compute_hessU_vec_prod,
       get_momentum_sample = function() L %*% rnorm(D))

}
