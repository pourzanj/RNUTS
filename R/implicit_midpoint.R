#' Single step implicit midpoint.
#'
#' Given a starting point, takes a single step using implicit midpoint.
#'
#' This function requires a function argument called do_nonlin_solve to solve
#' the nonlinear system that is required of implicit midpoint. do_nonlin_solve
#' can be as simple as a vanilla Newton iteration or more elaborate
#'
#' @param z0
#' @param z_1
#' @param direction
#' @param ham_system
#'
#' @return
#' @export
#'
#' @examples
take_one_step_im <- function(z0, z_1 = NULL, z_2 = NULL, direction, ham_system, do_nonlin_solve) {

  # unpack z0
  q0 <- z0$q
  p0 <- z0$p
  h <- z0$h*direction
  D <- length(q0) # dimension of q

  # unpack functions from ham_system
  M <- ham_system$M
  GradU <- function(q) ham_system$compute_gradU(q)
  HessU <- function(q) ham_system$compute_hessU(q)
  HessU_vec_prod <- function(q, v) ham_system$compute_hessU_vec_prod(q, v)

  # setup nonlinear solve functions required for backward euler half step
  # this includes function, its jacobian funciton, and a jacobian-vector product function
  g <- function(q_half) q_half - q0 - (h/2)*solve(M, p0 - (h/2)*GradU(q_half))
  Jg <- function(q_half) diag(D) + (h^2/4)*solve(M, HessU(q_half))
  Jg_v <- function(q_half, v) v + (h^2/4)*solve(M,HessU_vec_prod(q_half,v))

  # if we have two prev soln points use it to form a better Newton init guess
  x0 <- q0
  if (!is.null(z_1)) {
    q_1 <- z_1$q
    x0 <- ifelse(abs(q0+q_1)/abs(q0) < 0.1, 0.0, q0)
  }

  # do solve and extract info about num grad evals taken, whether newtons converged, etc
  nonlin_solve_soln <- do_nonlin_solve(g, Jg, Jg_v, x0)

  num_hess_evals <- nonlin_solve_soln$num_hess_evals
  num_hess_vec_prod_evals <- nonlin_solve_soln$num_hess_vec_prod_evals
  num_newton_iters <- nonlin_solve_soln$num_newton_iters
  error <- ifelse(nonlin_solve_soln$did_converge, NA, "Newton")

  # use solve to set backward euler half step
  q_half <- nonlin_solve_soln$x
  p_half <- p0 - (h/2)*GradU(q_half)

  # do forward euler half step
  q1 <- q.half + (h/2)*solve(M, p_half)
  p1 <- p.half - (h/2)*GradU(q_half)

  # pack back in to z1 and return
  z1 <- list(q = q1, p = p1, h = z0$h)

  list(z1 = z1, error = error, num_newton_iters = num_newton_iters,
       num_hess_evals = num_hess_evals, num_hess_vec_prod_evals = num_hess_vec_prod_evals)
}

impmid <- function(q0, p0, M, pos.first = TRUE, N, h) {

  D <- length(q0)

  q <- matrix(NA, nrow = N+1, ncol = D)
  p <- matrix(NA, nrow = N+1, ncol = D)
  omega <- rep(NA, N+1)
  num.newton.iters <- rep(NA, N+1)
  num_grad <- rep(NA,N+1)

  q[1,] <- q0
  p[1,] <- p0
  omega[1] <- HessU(q[1,]) %>% eigen %>% .$values %>% max %>% sqrt

  H0 <- (1/2)*sum(p0^2) + U(q0)
  for(n in 1:N) {

    if(pos.first) {

      # backward euler half step is reduced to just a nonlinear solve for p.half
      g <- function(q.half) q.half - q[n,] - (h/2)*solve(M,p[n,]) + (h^2/4)*solve(M,GradU(q.half))
      J <- function(q.half) diag(D) + (h^2/4)*solve(M,HessU(q.half))
      #newton.iter.soln <- newton_iteration(g, J, x0 = c(p[n,1],0.0), TOL = 1e-15,  max_iter = 100)
      x0 <- q[n,]
      # if(n >= 2) {
      #   #browser()
      #   x0 <- ifelse(abs(q[n,]+q[n-1,])/abs(q[n,]) < 0.1, 0.0, q[n,])
      # }
      newton.iter.soln <- newton_iteration(g, J, x0 = x0, TOL = 1e-13,  max_iter = 1000)
      #newton.iter.soln <- nleqslv(q[n,],g)
      num.newton.iters[n+1] <- newton.iter.soln$num.newton.iters
      num_grad[n+1] <- newton.iter.soln$num_grad
      #num.newton.iters[n+1] <- newton.iter.soln$nfcnt + D*newton.iter.soln$njcnt
      q.half <- newton.iter.soln$x1
      #q.half <- newton.iter.soln$x
      p.half <- p[n,] - (h/2)*GradU(q.half)

    } else {
      # backward euler half step is reduced to just a nonlinear solve for p.half
      g <- function(p.half) p.half - p[n,] + (h/2)*GradU(q[n,] + (h/2)*solve(M,p.half))
      J <- function(p.half) diag(D) + (h^2/4)*solve(M,HessU(q[n,] + (h/2)*solve(M,p.half)))
      #newton.iter.soln <- newton_iteration(g, J, x0 = c(p[n,1],0.0), TOL = 1e-15,  max_iter = 100)
      newton.iter.soln <- newton_iteration(g, J, x0 = p[n,], TOL = 1e-13,  max_iter = 100)
      num.newton.iters[n+1] <- newton.iter.soln$num.newton.iters
      p.half <- newton.iter.soln$x1
      q.half <- q[n,] + (h/2)*solve(M,p.half)
    }


    # forward euler half step
    q[n+1,] <- q.half + (h/2)*solve(M,p.half)
    p[n+1,] <- p.half - (h/2)*GradU(q.half)

    omega[n+1] <- HessU(q[n+1,]) %>% eigen %>% .$values %>% max %>% sqrt
  }

  H0 <- (1/2)*t(p0) %*% solve(M,p0) + U(q0)
  H0 <- as.vector(H0)
  H <- (1/2)*diag(p %*% solve(M,t(p))) + apply(q,1,U)
  H <- as.vector(H)
  soln <- tibble(t = c(0,cumsum(rep(h,N)))) %>%
    bind_cols(as_tibble(q)) %>%
    bind_cols(as_tibble(p)) %>%
    set_names(c("t", paste0("q",1:D), paste0("p",1:D))) %>%
    mutate(omega = omega) %>%
    mutate(num.newton.iters = num.newton.iters) %>%
    mutate(num_grad = num_grad) %>%
    mutate(H = H) %>%
    mutate(DeltaH = H0 - H) %>%
    mutate(accept.prob = pmin(1,exp(DeltaH)))

  return(soln)
}

