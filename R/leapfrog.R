take_one_step_lf <- function(z0, z_1, z_2, direction, ham_system, H0) {

  # unpack z0 and Hamiltonian system
  q0 <- z0$q
  p0 <- z0$p
  h <- z0$h*direction

  GradU <- function(q) ham_system$compute_gradU(q)
  M <- ham_system$M

  # leapfrog step
  p_half <- p0 - (h/2)*GradU(q0)
  q1 <- q0 + h*solve(M, p_half)
  p1 <- p_half - (h/2)*GradU(q1)

  # pack back in to z1
  z1 <- list(q = q1, p = p1, h = z0$h)

  # check for divergence
  H1 <- ham_system$compute_H(z1)
  error <- as.character(NA)
  if (H1-H0 >= 10) {
    error <- "Divergence"
    print("Divergence")
  }

  # return
  list(z1 = z1, error = error, num_grad = 2)
}
