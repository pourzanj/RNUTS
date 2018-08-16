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

get_nuts_samples <- function(num_samples, q0, h0, ham_system, integrator, DEBUG = FALSE) {

  D <- length(q0)
  q <- matrix(NA, nrow = num_samples+1, ncol = D)

  q[1,] <- q0

  lf <- create_integrator(is_implicit = FALSE)
  im <- create_integrator(is_implicit = TRUE)

  tree_depth <- rep(0, num_samples+1)
  total_grad_evals <- rep(0, num_samples+1)
  no_error <- rep(TRUE, num_samples+1)
  hist <- list(NA)

  for(iter in 1:num_samples){


    sample <- get_single_nuts_sample(as.vector(q[iter,]), p0 = NULL, h0, ham_system, integrator, max_treedepth = 10, DEBUG)
    q[iter+1,] <- sample$q

    tree_depth[iter+1] <- sample$hist %>% pull(depth) %>% max(na.rm = TRUE)
    total_grad_evals[iter+1] <- sample$hist %>% pull(num_grad) %>% sum(na.rm = TRUE)
    no_error[iter+1] <- (sample$hist$invalid == "U-Turn") %>% all(na.rm = TRUE)
    hist <- c(hist, list(sample$hist))
  }

  as_tibble(q) %>%
    set_names(paste0("q",1:D)) %>%
    mutate(tree_depth = tree_depth,
           total_grad_evals = total_grad_evals,
           no_error = no_error) %>%
    mutate(hist = hist)

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
    p0 <- as.vector(ham_system$get_momentum_sample())
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

  return(list(q = tree$z_rep$q, hist = tree$hist))
}