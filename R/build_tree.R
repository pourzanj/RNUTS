#' Create a history tibble for a onenode tree.
#'
#' A helper function to create_onenode_tree to create clean,
#' uniform tibbles that keep track of the history of the tree.
#'
#' @param z a list containing vectors q and p and a stepsize h representing a point in phase space plus a stepsize.
#' @param depth the depth in the whole NUTS tree that this node sits at
#' @param H the Hamiltonian value at this point in phase space
#' @param valid_subtree whether the subtree at that depth was valid
#' @param uturn whether a uturn occurred at that node
#' @param integrator_error either NA, "divergence", or "newton"
#' @param num_grad number of likelihood gradient evaluations it took to get that step
#' @param num_hess number of likelihood hessian evaluations it took to get that step
#' @param num_hess_vec number of likelihood hessian-vector product evaluations it took to get that step
#' @param num_newton  number of Newton iterations it took to get that step
#'
#' @return A tibble with a single row that represents a node in the tree including its depth, energy value, position, and whether the step was invalid
#' @export
#'
#' @examples
create_onenode_hist <- function(z, depth, H0, H, valid_subtree, uturn, integrator_error, num_grad, num_hess, num_hess_vec, num_newton) {
  D <- length(z$q)
  q <- matrix(z$q, nrow = 1) %>% as_tibble %>% set_names(paste0("q",1:D))
  p <- matrix(z$p, nrow = 1) %>% as_tibble %>% set_names(paste0("p",1:D))

  tibble(depth = depth, h = z$h, H0 = H0, H = H,
         valid_subtree = valid_subtree,
         uturn = uturn,
         integrator_error = integrator_error) %>%
    mutate(num_grad = num_grad,
           num_hess = num_hess,
           num_hess_vec = num_hess_vec,
           num_newton = num_newton) %>%
    bind_cols(bind_cols(q, p))
}

#' Create one node tree.
#'
#' This function is akin to a constructor. It makes sure
#' we create trees that uniformly have the same entries with the same names.
#'
#' @param z a list containing vectors q and p and a stepsize h representing a point in phase space plus a stepsize.
#' @param depth the depth in the whole NUTS tree that this node sits at
#' @param H the Hamiltonian value at this point in phase space
#' @param valid_subtree whether the subtree at that depth was valid
#' @param uturn whether a uturn occurred at that node
#' @param integrator_error either NA, "divergence", or "newton"
#' @param num_grad number of likelihood gradient evaluations it took to get that step
#' @param num_hess number of likelihood hessian evaluations it took to get that step
#' @param num_hess_vec number of likelihood hessian-vector product evaluations it took to get that step
#' @param num_newton  number of Newton iterations it took to get that step
#' @param DEBUG if this is on a tibble keeping track of the history of the node will be returned as well
#'
#' @return A one node tree which is eseentially a list which several attributes such as depth and whether the tree is valid.
#' @export
#'
#' @examples
create_onenode_tree <- function(z, depth, H0, H, valid_subtree, uturn, integrator_error, num_grad, num_hess, num_hess_vec, num_newton, DEBUG) {
  hist <- NULL
  if (DEBUG) {
    hist <- create_onenode_hist(z, depth, H0, H, valid_subtree, uturn, integrator_error, num_grad, num_hess, num_hess_vec, num_newton)
  }

  list(depth = depth,
       valid = valid_subtree,
       integrator_error = integrator_error,
       coordinate_uturn = rep(FALSE, length(z$q)),
       log_w = H0-H,
       rho = z$p,
       z_rep = z,
       z_minus = z,
       z_minus_1 = NULL,
       z_minus_2 = NULL,
       z_plus = z,
       z_plus_1 = NULL,
       z_plus_2 = NULL,
       num_grad = num_grad,
       num_hess = num_hess,
       num_hess_vec = num_hess_vec,
       num_newton = num_newton,
       hist = hist)
}

#' Build tree
#'
#' Build a NUTS tree starting from z0. If depth is 0, then this is just a single node.
#' If depth is 1, then it's two nodes. Depth is j then 2^j nodes. Tree are built recursively
#' e.g. if we need a depth 2 tree which has 4 nodes we'll build 2 trees of depth 1 and join
#' them together.
#'
#' Sometimes a tree can be invalid, either because there was a problem with the
#' integrator or because a U-Turn was detected. In this case the tree is marked as invalid and building
#' of the tree ceases.
#'
#' @param z0 Initial point to start from. Should contain q, p, h
#' @param z_1 z_{-1} Previous point. Useful for determining a guess of z_{1}
#' @param z_1 z_{-2} Previous, previous point. Useful for determining a guess of z_{1}
#' @param depth Number of levels of tree.
#' @param direction Direction we'd like to build tree in (forwards or backwards)
#' @param integrate_step a function that integrates a single step
#' @param DEBUG Flag to determine whether we return tibble that includes history of points
#'
#' @return A list reprenting a tree
#' @export
#'
#' @examples
build_tree <- function(depth, z0, z_1, z_2, direction, ham_system, H0, integrate_step, DEBUG = FALSE) {
  new_tree <- NULL

  # base case (take a single step)
  if(depth == 0){
    integrator_result <- integrate_step(z0, z_1, z_2, direction, ham_system, H0)
    new_tree <- create_onenode_tree(z = integrator_result$z1,
                                    depth = depth,
                                    H0 = H0,
                                    H = ham_system$compute_H(integrator_result$z1),
                                    valid_subtree = is.na(integrator_result$integrator_error),
                                    uturn = FALSE,
                                    integrator_error = integrator_result$integrator_error,
                                    num_grad = integrator_result$num_grad,
                                    num_hess = integrator_result$num_hess,
                                    num_hess_vec = integrator_result$num_hess_vec,
                                    num_newton = integrator_result$num_newton,
                                    DEBUG = DEBUG)
  }

  # recursion
  else{

    inner_subtree <- build_tree(depth-1, z0, z_1, z_2, direction, ham_system, H0, integrate_step, DEBUG)

    # only build outer subtree and tack it on if inner subtree was valid. otherwise
    # just return the inner_subtree
    if (inner_subtree$valid) {

      # assume direciton is forward in which case we build outer subtree starting
      # from z_plus. If direction is backward then we start from z_minus
      z0.outer <- inner_subtree$z_plus
      if (direction == -1) {
        z0.outer <- inner_subtree$z_minus
      }

      # build outer_subtree and join even if it's invalid because we might
      # want to view its history.
      outer_subtree <- build_tree(depth-1, z0.outer, z_1, z_2, direction, ham_system, H0, integrate_step, DEBUG)
      new_tree <- join_subtrees(inner_subtree, outer_subtree, direction, biased_progressive_sampling = FALSE, ham_system, DEBUG)

    } else {
      new_tree <- inner_subtree
    }

  }

  # if we're in debug mode update the depth of the subtree in the history tibble. note that this has to be
  # done here in build_tree because join_trees is also used in the main NUTS call where we don't want these depths to be updated in this manner
  if (DEBUG) {
    depth_ <- depth
    new_tree$hist <- new_tree$hist %>%
      mutate(depth = depth_) %>%
      mutate(valid_subtree = new_tree$valid)
  }

  new_tree
}


#' Join two subtrees.
#'
#' Joins an inner subtree and an outer subtree.
#'
#' This function is called either in sampling where we tack on a new subtree to our current tree
#' The other place where it's called is in build_tree() where we create two trees and join them.
#'
#' We can always assume both trees passed in are non-null but we can't always assume they're the
#' same size. We can assume the inner_tree is valid, but not the outer. The outer subtree could have stopped early because it was invalid. In which case
#' the representative sample will just be from the inner_subtree
#'
#' The represenative sample is in in accordance with Betancourt's continuous sampling of subtrees,
#' where every tree has a representative node z_rep. For a 1-node tree, the representative is that one node.
#' For a two node tree we sample a representative from the two nodes. This sampling is
#' described in the function sample_new_represenative().
#'
#' @param inner_subtree The valid and non-null subtree created first.
#' @param outer_subtree The subtree created second as an extension of the end of the inner subtree.
#' @param direction the direction the tree were build determines if outer is added on the plus side or minus side
#' @param ham_system the hamiltonian system of the problem
#' @param biased_progressive_sampling whether to sample biasedly sample the new tree (should only be done in the original call to build_tree() called from NUTS)
#' @param DEBUG whether to include the history tibble
#'
#' @return A new joined tree.
#' @export
#'
#' @examples
join_subtrees <- function(inner_subtree, outer_subtree, direction, biased_progressive_sampling, ham_system, DEBUG) {

  # returned tree starts as copy of inner_subtree.
  tree <- inner_subtree
  tree$depth <- outer_subtree$depth + 1
  tree$log_w <- log(exp(inner_subtree$log_w) + exp(outer_subtree$log_w))

  # only update representative if outer_subtree is valid
  if(outer_subtree$valid) {
    if(biased_progressive_sampling) {
      tree$z_rep <- sample_new_representative_biasedly(inner_subtree, outer_subtree)
    } else{
      tree$z_rep <- sample_new_representative_uniformly(inner_subtree, outer_subtree)
    }

  }

  # update z_plus, z_plus_1, z_plus_2, z_minus, z_minus_1, z_minus_2. note if depth of the
  # new joined tree == 1 then these are set to NULL so we manually set them. if depth == 2
  # we'll have four nodes but neither of the subtree will join will have z_plus_2 and z_minus_2
  # so we need to set those manually as well
  if (direction == 1) {
    tree$z_plus = outer_subtree$z_plus
    tree$z_plus_1 = outer_subtree$z_plus_1
    tree$z_plus_2 = outer_subtree$z_plus_2
  } else {
    tree$z_minus = outer_subtree$z_minus
    tree$z_minus_1 = outer_subtree$z_minus_1
    tree$z_minus_2 = outer_subtree$z_minus_2
  }

  if (tree$depth == 1) {
    tree$z_plus_1 = tree$z_minus
    tree$z_minus_1 = tree$z_plus
  }

  if (tree$depth == 2) {
    if (direction == 1) {
      tree$z_plus_2 = inner_subtree$z_plus
      tree$z_minus_2 = outer_subtree$z_minus
    } else {
      tree$z_plus_2 = outer_subtree$z_plus
      tree$z_minus_2 = inner_subtree$z_minus
    }
  }

  # check to see if new joined tree is valid. This only happens if both subtrees
  # are valid and if the u-turn criteria is met between z_plus and z_minus
  both_subtrees_valid <- inner_subtree$valid & outer_subtree$valid
  tree$coordinate_uturn <- update_coordinate_uturn(tree, ham_system)
  nouturn_criteria_met <- check_uturn_criteria(tree, ham_system)
  #nouturn_criteria_met <- !all(tree$coordinate_uturn)
  if(both_subtrees_valid & nouturn_criteria_met) {
    tree$valid = TRUE
  } else {
    tree$valid = FALSE
  }

  # update the integrator error if there was one
  if(!outer_subtree$valid & !is.na(outer_subtree$integrator_error)) {
    tree$integrator_error <- outer_subtree$integrator_error
  } else {
    tree$integrator_error <- as.character(NA)
  }

  # update number of evals
  tree$num_grad <- sum(inner_subtree$num_grad, outer_subtree$num_grad, na.rm = TRUE)
  tree$num_hess <- sum(inner_subtree$num_hess, outer_subtree$num_hess, na.rm = TRUE)
  tree$num_hess_vec <- sum(inner_subtree$num_hess_vec, outer_subtree$num_hess_vec, na.rm = TRUE)
  tree$num_newton <- sum(inner_subtree$num_newton, outer_subtree$num_newton, na.rm = TRUE)

  # update hist if we're in debug mode
  if (DEBUG) {
    tree$hist <- join_tree_histories(inner_subtree, outer_subtree, direction, nouturn_criteria_met, both_subtrees_valid)
  }

  tree
}

join_tree_histories <- function(inner_subtree, outer_subtree, direction, nouturn_criteria_met, both_subtrees_valid) {

  new_hist <- NULL

  # append histograms in the order that depends on which direciton we're going

  # if the noturn criteria was NOT met AND both subtrees were valid then that means
  # this was the subtree that caused the U-Turn indicator to go off. If one of the trees
  # was already invalid then this U-Turn wasn't what causes the tree to go invalid
  # because it already was before due to an earlier subtree
  if (direction == 1) {
    new_hist <- bind_rows(inner_subtree$hist, outer_subtree$hist)
  } else {
    new_hist <- bind_rows(outer_subtree$hist, inner_subtree$hist)
  }

  # if the noturn criteria was NOT met AND both subtrees were valid then that means
  # this was the subtree that caused the U-Turn indicator to go off. If one of the trees
  # was already invalid then this U-Turn wasn't what causes the tree to go invalid
  # because it already was before due to an earlier subtree
  if (!nouturn_criteria_met & both_subtrees_valid) {
    new_hist <- new_hist %>% mutate(valid_subtree = FALSE)
    new_hist$uturn[1] <- TRUE
    new_hist$uturn[nrow(new_hist)] <- TRUE
  }

  new_hist
}

#' Take a uniform sample over the joined trajectory.
#'
#' Uniformly sample the joined trajectory according the weights of the two constituent parts.
#'
#' @param inner_subtree a valid inner_subtree
#' @param outer_subtree a valid outer_Subtree
#'
#' @return either representative sample from the left or right subtree
#' @export
#'
#' @examples
sample_new_representative_uniformly <- function(inner_subtree, outer_subtree) {

  new_z_rep <- inner_subtree$z_rep

  log_w_old <- inner_subtree$log_w
  log_w_new <- outer_subtree$log_w
  log_w_total <- log(exp(log_w_old) + exp(log_w_new))

  # if the new tree has greater weight then sample it with prob. 1
  # otherwise sample with a prob. proportional the respect tree weights
  if (log_w_new >= log_w_total) {
    new_z_rep <- outer_subtree$z_rep
  } else {
    # sample a bernoulli with prob. that depends on weight
    # if TRUE then we take rep from outer_subtree else the rep remains the
    # the rep from the old (inner) subtree
    if (runif(1) <= exp(log_w_new - log_w_total)) {
      new_z_rep <- outer_subtree$z_rep
    }
  }

  new_z_rep
}

#' Take a uniform sample over the joined trajectory.
#'
#' If the outer tree has a bigger weight then we use it's representative sample
#' with probability one.
#' Otherwise, we'll make the new representative either the representative of the old
#' tree with probability proportional to the weight of the old tree and we'll make the
#' new representative the representative of the new subtree with probability proportional
#' to the weight of the new subtree.
#'
#' @param inner_subtree a valid inner_subtree
#' @param outer_subtree a valid outer_Subtree
#'
#' @return either representative sample from the left or right subtree
#' @export
#'
#' @examples
sample_new_representative_biasedly <- function(inner_subtree, outer_subtree) {

  new_z_rep <- inner_subtree$z_rep

  log_w_old <- inner_subtree$log_w
  log_w_new <- outer_subtree$log_w

  # if the new tree has greater weight then sample it with prob. 1
  # otherwise sample with a prob. proportional the respect tree weights
  if (log_w_new >= log_w_old) {
    new_z_rep <- outer_subtree$z_rep
  } else {
    # sample a bernoulli with prob. that depends on weight
    # if TRUE then we take rep from outer_subtree else the rep remains the
    # the rep from the old (inner) subtree
    if (runif(1) <= exp(log_w_new - log_w_old)) {
      new_z_rep <- outer_subtree$z_rep
    }
  }

  new_z_rep
}

#' Check that there's NO U-Turns
#'
#' Returns true if there's no U-Turn in either direction.
#'
#' @param tree the joined tree we're checking the U-Turn of
#' @param ham_system the Hamiltonian system for the problem
#'
#' @return TRUE if there's no U-Turn and FALSE if there is
#' @export
#'
#' @examples
check_uturn_criteria <- function(tree, ham_system) {

  q_plus <- tree$z_plus$q
  q_minus <- tree$z_minus$q

  # instead of momentums get velocities by multipling by M^{-1}
  M <- ham_system$M
  # v_plus <- solve(M,tree$z_plus$p)
  # v_minus <- solve(M,tree$z_minus$p)
  v_plus <- tree$z_plus$p
  v_minus <- tree$z_minus$p

  no_uturn_forward <- as.numeric(v_plus %*% (q_plus-q_minus)) > 0
  no_uturn_backward <- as.numeric(-v_minus %*% (q_minus-q_plus)) > 0

  no_uturn_forward & no_uturn_backward
}

#' Check that there's NO U-Turns using new generalized Criteria.
#'
#' Returns true if there's no U-Turn in either direction.
#'
#' @param tree the joined tree we're checking the U-Turn of
#' @param ham_system the Hamiltonian system for the problem
#'
#' @return TRUE if there's no U-Turn and FALSE if there is
#' @export
#'
#' @examples
check_generalized_uturn_criteria <- function(tree, ham_system) {

  q_plus <- tree$z_plus$q
  q_minus <- tree$z_minus$q

  # instead of momentums get velocities by multipling by M^{-1}
  M <- ham_system$M
  v_plus <- tree$z_plus$p
  v_minus <- tree$z_minus$p

  no_uturn_forward <- as.numeric(v_plus %*% (q_plus-q_minus)) > 0
  no_uturn_backward <- as.numeric(-v_minus %*% (q_minus-q_plus)) > 0

  no_uturn_forward & no_uturn_backward
}

#' Check for U-Turns at the coordinate level
#'
#' In each coordinate q(t)-q(0) tells us whether we've gone in the negative or
#' or positive direction. If we multiuply that by p(t) for that dimension we get
#' whether we're U-Turning specifically in that dimension. We keep track of this because
#' it may be a good U-Turn criteria to check that ALL dimensions have U-Turned
#'
#' @param tree the newly joined tree we're going to check
#' @param ham_system the Hamiltonian system for the problem
#'
#' @return a vector tracking the U-Turn status of each coordinate
#' @export
#'
#' @examples
update_coordinate_uturn <- function(tree, ham_system) {

  q_plus <- tree$z_plus$q
  q_minus <- tree$z_minus$q

  # instead of momentums get velocities by multipling by M^{-1}
  M <- ham_system$M
  v_plus <- solve(M,tree$z_plus$p)
  v_minus <- solve(M,tree$z_minus$p)

  # this is a vector that is true if that dimension U-Turned
  coordinate_uturns_forward <- (v_plus*(q_plus-q_minus) < 0)
  coordinate_uturns_backward <- (-v_minus*(q_minus-q_plus) < 0)

  tree$coordinate_uturn | coordinate_uturns_forward | coordinate_uturns_backward
}
