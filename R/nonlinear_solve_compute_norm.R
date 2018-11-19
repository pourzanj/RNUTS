compute_two_norm <- function(v) {
  sqrt(sum(v^2))
}

compute_rms_norm <- function(v) {
  sqrt(mean(v^2))
}

create_compute_weighted_rms_norm <- function(rtols, atols) {
  function(v) {
    # weights
    w <- rtols*abs(v) + atol
    sqrt(mean((v/w)^2))
  }
}
