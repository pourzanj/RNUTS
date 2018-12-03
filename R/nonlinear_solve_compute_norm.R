compute_two_norm <- function(v) {
  sqrt(sum(v^2))
}

compute_rms_norm <- function(v) {
  sqrt(mean(v^2))
}

#' Create a function to compute Weighted RMS Norm.
#'
#' @param w approximate scale of each variable
#'
#' @return a function that computes a weighted RMS norm on a vector
#' @export
#'
#' @examples
create_compute_weighted_rms_norm <- function(w) {
  function(v) {
    # weights
    sqrt(mean((v/w)^2))
  }
}
