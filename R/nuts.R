create_hamiltonian_system <- function(M, compute_U, compute_gradU) {

  if (is.matrix(M)) {
    stop("Only diaganol mass matrices specified as vectors are supported as of now.")
  }

  D <- length(M)
  L <- chol(diag(M))

  list(M = M,
       compute_gradU = compute_gradU,
       compute_H = function(z) (1/2)*sum(z$p*(1/M)*z$p) + compute_U(z$q),
       get_momentum_sample = function() L %*% rnorm(D))

}

get_nuts_samples <- function(num_samples, U, GradU, q0, h0, Sigma, implicit = FALSE, uturn_criteria, min_DeltaH, DEBUG = FALSE) {

  D <- length(q0)
  q <- matrix(NA, nrow = num_samples+1, ncol = D)

  chol_Sigma <- chol(Sigma)

  q[1,] <- q0

  if(DEBUG) {

    tree.depth <- rep(0, num_samples+1)
    total.grad.evals <- rep(0, num_samples+1)
    accepted <- rep(TRUE, num_samples+1)
    divergent <- rep(FALSE, num_samples+1)
    newton.failure <- rep(FALSE, num_samples+1)
    hist <- list(NA)

    for(iter in 1:num_samples){
      nuts <- NUTS_one_step(U, GradU, q[iter,], p0 = NULL, h0, chol_Sigma, implicit, uturn_criteria, min_DeltaH, max_treedepth = 10, DEBUG = TRUE)
      q[iter+1,] <- nuts$q
      tree.depth[iter+1] <- nuts$hist %>% pull(depth) %>% max(na.rm = TRUE)
      total.grad.evals[iter+1] <- nuts$hist %>% pull(NumGradEval) %>% sum(na.rm = TRUE)
      accepted[iter+1] <- ifelse(sum(as.vector(nuts$q) == q[iter,]), FALSE, TRUE)
      divergent[iter+1] <- (nuts$hist %>% pull(divergent) %>% sum(na.rm = TRUE) >= 1)
      newton.failure[iter+1] <- (nuts$hist %>% pull(newton.failure) %>% sum(na.rm = TRUE) >= 1)
      hist <- c(hist, list(nuts$hist))
    }

    return(samples = as_tibble(q) %>% set_names(paste0("q",1:D)) %>%
             mutate(tree.depth = tree.depth, total.grad.evals = total.grad.evals, accepted = accepted,
                    divergent = divergent, newton.failure = newton.failure) %>%
             mutate(hist =hist))

  } else {
    for(iter in 1:num_samples){

      nuts <- NUTS_one_step(U, GradU, q[iter,], p0 = NULL, h0, chol_Sigma, implicit, uturn_criteria, min_DeltaH, max_treedepth = 10, DEBUG = FALSE)
      q[iter+1,] <- nuts$q

    }

    return(as_tibble(q) %>% set_names(paste0("q",1:D)))
  }


}

#' Get single NUTS sample
#'
#' @param z0
#' @param integrator
#' @param max_treedepth
#' @param DEBUG
#'
#' @return List containing sample and if debugging a tibble representing the history
#' @export
#'
#' @examples
get_single_nuts_sample <- function(q0, p0, h0, ham_system, integrator, max_treedepth = 10, DEBUG = FALSE) {

  # sample momentum
  if (is.null(p0)) {
    p0 <- ham_system$get_momentum_sample()
  }
  z0 <- list(q = q0, p = p0, h = h0)

  # start building tree
  tree <- create_onenode_tree(depth = NA,
                              invalid = NA,
                              z = z0,
                              ham_system = ham_system,
                              DEBUG = DEBUG)

  # sample directions we'll go in ahead of time for easier debugging
  directions <- base::sample(c(-1, 1), max_treedepth, replace = TRUE)
  for(depth in 0:(max_treedepth-1)) {

    new_subtree <- NULL

    # we can either evolve the right-most node right (z_plus), or the left-most node left (z_minus)
    if(directions[depth+1] == 1){
      new_subtree <- build_tree(depth, tree$z_plus, tree$z_plus_1, directions[depth+1], ham_system, integrator, DEBUG)
    }
    else{
      new_subtree <- build_tree(depth, tree$z_minus, tree$z_minus_1, directions[depth+1], ham_system, integrator, DEBUG)
    }

    tree <- join_subtrees(tree, new_subtree, directions[depth+1], DEBUG)
    if (tree$invalid) break
  }

  return(list(q = tree$z_rep, hist = tree$hist))
}
