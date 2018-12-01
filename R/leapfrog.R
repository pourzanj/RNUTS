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
  if (H1-H0 >= 1000) {
    error <- "Divergence"
    warning("Divergence detected")
  }

  # return
  list(z1 = z1, error = error, num_grad = 2)
}

# leapfrog <- function(q0, p0, M, N, h) {
#
#   D <- length(q0)
#
#   q <- matrix(NA, nrow = N+1, ncol = D)
#   p <- matrix(NA, nrow = N+1, ncol = D)
#   omega <- rep(NA, N+1)
#
#   q[1,] <- q0
#   p[1,] <- p0
#   omega[1] <- HessU(q[1,]) %>% eigen %>% .$values %>% max %>% sqrt
#
#   # start solving
#   for(n in 1:N) {
#
#     phalf <- p[n,] - (h/2)*GradU(q[n,])
#     q[n+1,] <- q[n,] + h*solve(M,phalf)
#     p[n+1,] <- phalf - (h/2)*GradU(q[n+1,])
#     omega[n+1] <- HessU(q[n+1,]) %>% eigen %>% .$values %>% max %>% sqrt
#
#   }
#
#   H0 <- (1/2)*t(p0) %*% solve(M,p0) + U(q0)
#   H0 <- as.vector(H0)
#   H <- (1/2)*diag(p %*% solve(M,t(p))) + apply(q,1,U)
#   H <- as.vector(H)
#   soln <- tibble(t = c(0,cumsum(rep(h,N)))) %>%
#     bind_cols(as_tibble(q)) %>%
#     bind_cols(as_tibble(p)) %>%
#     set_names(c("t", paste0("q",1:D), paste0("p",1:D))) %>%
#     mutate(omega = omega) %>%
#     mutate(H = H) %>%
#     mutate(DeltaH = H0 - H) %>%
#     mutate(acceptprob = pmin(1,exp(DeltaH)))
#
#   return(soln)
# }
