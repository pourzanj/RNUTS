create_integrator <- function(is_implicit) {

  # set integrator to either explicit or implicit
  integrate_one_step <- take_one_step_lf
  if (is_implicit) {
    integrate_one_step <- take_one_step_im
  }

  # define and return step function
  take_step <- function(z0, z_1, direction, ham_system, DEBUG) {


    result <- integrate_one_step(z0, z_1, direction, ham_system)
    H <- ham_system$compute_H(result$z1)

    create_onenode_tree(depth = 0,
                        invalid = result$error,
                        z = result$z1,
                        ham_system = ham_system,
                        DEBUG = DEBUG)
  }

  list(take_step = take_step)
}



take_one_step_lf <- function(z0, z_1, direction, ham_system) {

  # unpack z0
  q0 <- z0$q
  p0 <- z0$p
  h <- z0$h*direction
  GradU <- function(q) ham_system$compute_gradU(q)

  # leapfrog step
  p_half <- p0 - (h/2)*GradU(q0)
  q1 <- q0 + h*(1/M)*p_half
  p1 <- p_half - (h/2)*GradU(q1)

  # pack back in to z1
  z1 <- list(q = q1, p = p1, h = z0$h)

  # check for divergence
  H0 <- ham_system$compute_H(z0)
  H1 <- ham_system$compute_H(z1)
  error <- ifelse(H0-H1 <= 1000, NA, "Divergence")

  # return
  list(z1 = z1, error = error)
}
