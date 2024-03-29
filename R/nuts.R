#' Get multiple NUTS samples from a posterior
#'
#' @param num_samples
#' @param q0
#' @param h0
#' @param ham_system
#' @param integrator
#' @param max_treedepth
#' @param DEBUG
#'
#' @return
#' @export
#'
#' @examples
get_nuts_samples <- function(num_samples, q0, h0, ham_system, integrator, refresh = NULL, max_treedepth = 10, DEBUG = FALSE) {

  if(is.null(refresh)) {
    refresh <- floor(num_samples/10)
  }

  D <- length(q0)
  q <- matrix(NA, nrow = num_samples+1, ncol = D)

  q[1,] <- q0

  tree_depth <- rep(as.integer(NA), num_samples+1)
  error <- rep(as.character(NA), num_samples+1)

  num_grad <- rep(as.numeric(NA), num_samples+1)
  num_hess <- rep(as.numeric(NA), num_samples+1)
  num_hess_vec <- rep(as.numeric(NA), num_samples+1)
  num_newton <- rep(as.numeric(NA), num_samples+1)

  hist <- list(NA)

  for(iter in 1:num_samples){

    sample <- get_single_nuts_sample(as.vector(q[iter,]), p0 = NULL, h0, ham_system, integrator, max_treedepth, DEBUG)
    q[iter+1,] <- sample$q

    tree_depth[iter+1] <- sample$tree_depth
    error[iter+1] <- sample$error

    num_grad[iter+1] <- sample$num_grad
    num_hess[iter+1] <- sample$num_hess
    num_hess_vec[iter+1] <- sample$num_hess_vec
    num_newton[iter+1] <- sample$num_newton

    if(DEBUG) {
      hist <- c(hist, list(sample$hist))
    }

    if(iter %% refresh == 0) {
      print(paste("Sample", iter, "completed. Average Depth: ", mean(tree_depth[(iter-refresh+1):iter])))
    }

  }

  tibble(tree_depth = tree_depth,
         error = error,
         num_grad = num_grad,
         num_hess = num_hess,
         num_hess_vec = num_hess_vec,
         num_newton = num_newton,
         hist = hist) %>%
         bind_cols(as_tibble(q) %>% set_names(paste0("q",1:D)))
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
  H0 <- ham_system$compute_H(z0)

  # start building tree
  tree <- create_onenode_tree(z = z0,
                              depth = NA,
                              H0 = H0,
                              H = H0,
                              valid_subtree = TRUE,
                              uturn = FALSE,
                              integrator_error = as.character(NA),
                              num_grad = as.integer(NA),
                              num_hess = as.integer(NA),
                              num_hess_vec = as.integer(NA),
                              num_newton = as.integer(NA),
                              DEBUG = DEBUG)

  # sample directions we'll go in ahead of time for easier debugging
  directions <- base::sample(c(-1, 1), max_treedepth, replace = TRUE)
  for(depth in 0:(max_treedepth-1)) {
    new_subtree <- NULL

    # we can either evolve the right-most node right (z_plus), or the left-most node left (z_minus)
    if(directions[depth+1] == 1){
      new_subtree <- build_tree(depth, tree$z_plus, tree$z_plus_1, tree$z_plus_2, directions[depth+1], ham_system, H0, integrator, DEBUG)
    }
    else{
      new_subtree <- build_tree(depth, tree$z_minus, tree$z_minus_1, tree$z_minus_2, directions[depth+1], ham_system, H0, integrator, DEBUG)
    }
    tree <- join_subtrees(tree, new_subtree, directions[depth+1], biased_progressive_sampling = TRUE, ham_system, DEBUG)
    if (!tree$valid) break
  }

  return(list(q = tree$z_rep$q, tree_depth = tree$depth,  error = tree$integrator_error,
              num_grad = tree$num_grad,
              num_hess = tree$num_hess,
              num_hess_vec = tree$num_hess_vec,
              num_newton = tree$num_newton,
              hist = tree$hist))
}
