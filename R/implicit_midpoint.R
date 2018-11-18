#' Title
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
take_one_step_im <- function(z0, z_1, direction, ham_system, nonlinear_solve) {

  # unpack z0
  q0 <- z0$q
  p0 <- z0$p
  h <- z0$h*direction
  GradU <- function(q) ham_system$compute_gradU(q)
  M <- ham_system$M

  # backward euler half step
  g <- function(q.half) q.half - q0 - (h/2)*(1/M)*p0 + (h^2/4)*(1/M)*GradU(q.half)
  gHessVecProd <- function(q.half, v) v + (h^2/4)*(1/M)*HessVecProd(q.half,v)
  Jg <- function(q.half) diag(length(q0)) + (h^2/4)*(1/M)*HessU(q.half)

  #newton.iter.soln <- InexactNewtonWithBacktracking(diag(M), h0,g, gHessVecProd, q0, nTOL = 1e-15)
  x0 <- q0
  if (!is.null(z_1)) {
    q_1 <- z_1$q
    x0 <- ifelse(abs(q0+q_1)/abs(q0) < 0.1, 0.0, q0)
  }

  newton.iter.soln <- newton_iteration(g, Jg, x0 = x0, TOL = 1e-13, max_iter = 100)
  num_grad <- newton.iter.soln$num_grad
  num_newton <- newton.iter.soln$num.newton.iters
  newton.converged <- newton.iter.soln$converged

  q.half <- newton.iter.soln$x
  p.half <- p0 - (h/2)*GradU(q.half)

  # forward euler half step
  q1 <- q.half + (h/2)*(1/M)*p.half
  p1 <- p.half - (h/2)*GradU(q.half)

  # pack back in to z1
  z1 <- list(q = q1, p = p1, h = z0$h)

  # check for Newton not Converging
  error <- ifelse(newton.converged, NA, "Newton")

  # return
  list(z1 = z1, error = error, num_newton = num_newton, num_grad = num_grad)
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

