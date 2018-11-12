#' Title
#'
#' @param z
#' @param depth
#' @param H
#' @param invalid
#' @param num_newton
#' @param num_grad
#'
#' @return
#' @export
#'
#' @examples
create_onenode_hist <- function(z, depth, H, invalid, num_newton, num_grad) {

  if(!is.character(invalid)) {
    invalid <- as.character(NA)
  }
  D <- length(z$q)
  q <- matrix(z$q, nrow = 1) %>% as_tibble %>% set_names(paste0("q",1:D))
  p <- matrix(z$p, nrow = 1) %>% as_tibble %>% set_names(paste0("p",1:D))

  tibble(depth = depth, H = H, invalid = invalid) %>%
    mutate(num_newton = num_newton, num_grad = num_newton) %>%
    bind_cols(bind_cols(q, p))
}

#' Title
#'
#' use this function to create trees to make sure all trees have the same entries with the same names
#'
#' @param depth
#' @param invalid
#' @param z
#' @param ham_system
#' @param num_newton
#' @param num_grad
#' @param DEBUG
#'
#' @return
#' @export
#'
#' @examples
create_onenode_tree <- function(depth, invalid, z, ham_system, num_newton = NA, num_grad = NA, DEBUG) {

  H <- ham_system$compute_H(z)

  hist <- NULL
  if (DEBUG) {
    hist <- create_onenode_hist(z, depth, H, invalid, num_newton, num_grad)
  }

  list(depth = depth,
       invalid = !is.na(invalid),
       coordinate_uturn = rep(FALSE, length(z$q)),
       log_w = -H,
       z_rep = z,
       z_minus = z,
       z_minus_1 = NULL,
       z_plus = z,
       z_plus_1 = NULL,
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
#' @param depth Number of levels of tree.
#' @param direction Direction we'd like to build tree in (forwards or backwards)
#' @param integrator List representing our integrator
#' @param DEBUG Flag to determine whether we return tibble that includes history of points
#'
#' @return A list reprenting a tree
#' @export
#'
#' @examples
build_tree <- function(depth, z0, z_1, direction, ham_system, integrator, DEBUG = FALSE) {

  new_tree <- NULL

  # base case (take a single step)
  if(depth == 0){
    new_tree <- integrator$take_step(z0, z_1, direction, ham_system, DEBUG)
  }

  # recursion
  else{

    inner_subtree <- build_tree(depth-1, z0, z_1, direction, ham_system, integrator, DEBUG)

    # only build outer subtree and tack it on if inner subtree was valid. otherwise
    # just return the inner_subtree
    if (!inner_subtree$invalid) {

      # assume direciton is forward in which case we build outer subtree starting
      # from z_plus. If direction is backward then we start from z_minus
      z0.outer <- inner_subtree$z_plus
      if (direction == -1) {
        z0.outer <- inner_subtree$z_minus
      }

      outer_subtree <- build_tree(depth-1, z0.outer, z_1, direction, ham_system, integrator, DEBUG)
      new_tree <- join_subtrees(inner_subtree, outer_subtree, direction, DEBUG)
    } else {
      new_tree <- inner_subtree
    }

  }

  # depth has to be set manually when we're debugging
  if (DEBUG) {
    d <- depth
    new_tree$hist <- new_tree$hist %>% mutate(depth = d)
  }

  new_tree
}


#' Join Subtrees
#'
#' This function is called either in sampling where we tack on a new subtree to our current tree
#' The other place where it's called is in build_tree() where we create two trees and join them.
#'
#' We can always assume both trees passed in are non-null but we can't always assume they're the
#' same size. We also can't assume they're both valid. The outer subtree could have stopped early because it was invalid. In which case
#' the representative sample will just be from the inner_subtree
#'
#' The represenative sample is in in accordance with Betancourt's continuous sampling of subtrees,
#' where every tree has a representative node z_rep. For a 1-node tree, the representative is that one node.
#' For a two node tree we sample a representative from the two nodes. This sampling is
#' described in the function sample_new_represenative().
#'
#'
#' @param inner_subtree
#' @param outer_subtree
#'
#' @return
#' @export
#'
#' @examples
join_subtrees <- function(inner_subtree, outer_subtree, direction, DEBUG) {

  # returned tree starts as copy of inner_subtree. the depths of the two trees
  # being joined should always be the same and depth of returned tree should be +1
  tree <- inner_subtree
  tree$depth <- outer_subtree$depth + 1
  tree$log_w <- log(exp(inner_subtree$log_w) + exp(outer_subtree$log_w))
  tree$z_rep <- sample_new_representative(inner_subtree, outer_subtree)

  # update z_plus, z_plus_1, z_minus, z_minus_1. note if depth == 1 then these are set
  # to NULL so we manually set them
  if (direction == 1) {
    tree$z_plus = outer_subtree$z_plus
    tree$z_plus_1 = outer_subtree$z_plus_1
  } else {
    tree$z_minus = outer_subtree$z_minus
    tree$z_minus_1 = outer_subtree$z_minus_1
  }

  if (tree$depth == 1) {
    tree$z_plus_1 = tree$z_minus
    tree$z_minus_1 = tree$z_plus
  }

  # check to see if new joined tree is valid. This only happens if both subtrees
  # are valid and if the u-turn criteria is met between z_plus and z_minus
  both_subtrees_valid <- !inner_subtree$invalid & !outer_subtree$invalid
  tree$coordinate_uturn <- update_coordinate_uturn(tree, system)
  nouturn_criteria_met <- check_uturn_criteria(tree, system)
  #nouturn_criteria_met <- !all(tree$coordinate_uturn)

  if(both_subtrees_valid & nouturn_criteria_met) {
    tree$invalid = FALSE
  } else {
    tree$invalid = TRUE
  }

  # update hist if we're in debug mode
  if (DEBUG) {
    tree$hist <- join_tree_histories(inner_subtree, outer_subtree, direction, nouturn_criteria_met, both_subtrees_valid)
  }

  tree
}

join_tree_histories <- function(inner_subtree, outer_subtree, direction, nouturn_criteria_met, both_subtrees_valid) {

  new_hist <- NULL

  # append histograms in the order that depends on which direciton we're going
  if (direction == 1) {
    new_hist <- bind_rows(inner_subtree$hist, outer_subtree$hist)
  } else {
    new_hist <- bind_rows(outer_subtree$hist, inner_subtree$hist)
  }

  # if the noturn criteria was NOT met AND both subtrees were valid then that means
  # this was teh subtree that caused the U-Turn indicator to go off. If one of the trees
  # was already invalid then this U-Turn wasn't what causes the tree to go invalid
  # because it already was before due to an earlier subtree
  if (!nouturn_criteria_met & both_subtrees_valid) {
    new_hist <- new_hist %>% mutate(invalid = "U-Turn")
  }

  new_hist
}

#' Take biased tree sample
#'
#' First make sure outer subtree is valid. If it's not then the representative remains
#' the representatitive from the inner subtree.
#'
#' If the outer tree has a bigger weight then we use it's representative sample
#' with probability one.
#' Otherwise, we'll make the new representative either the representative of the old
#' tree with probability proportional to the weight of the old tree and we'll make the
#' new representative the representative of the new subtree with probability proportional
#' to the weight of the new subtree.
#'
#' @param inner_subtree
#' @param outer_subtree
#'
#' @return either representative sample from the left or right subtree
#' @export
#'
#' @examples
sample_new_representative <- function(inner_subtree, outer_subtree) {

  new_z_rep <- inner_subtree$z_rep

  if (!outer_subtree$invalid) {
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

  }

  new_z_rep
}

#' Check that there's NO U-Turns
#'
#' Returns true if there's no U-Turn in either direction
#'
#' @param tree
#' @param ham_system
#'
#' @return
#' @export
#'
#' @examples
check_uturn_criteria <- function(tree, ham_system) {

  q_plus <- tree$z_plus$q
  q_minus <- tree$z_minus$q

  # instead of momentums get velocities by multipling by M^{-1}
  v_plus <- (1/M)*tree$z_plus$p
  v_minus <- (1/M)*tree$z_minus$p

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
#' @param tree
#' @param ham_system
#'
#' @return
#' @export
#'
#' @examples
update_coordinate_uturn <- function(tree, ham_system) {

  q_plus <- tree$z_plus$q
  q_minus <- tree$z_minus$q

  # instead of momentums get velocities by multipling by M^{-1}
  v_plus <- (1/M)*tree$z_plus$p
  v_minus <- (1/M)*tree$z_minus$p

  # this is a vector that is true if that dimension U-Turned
  coordinate_uturns_forward <- (v_plus*(q_plus-q_minus) < 0)
  coordinate_uturns_backward <- (-v_minus*(q_minus-q_plus) < 0)

  tree$coordinate_uturn | coordinate_uturns_forward | coordinate_uturns_backward
}








