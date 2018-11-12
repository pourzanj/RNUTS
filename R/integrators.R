#' Title
#'
#' @param is_implicit
#'
#' @return
#' @export
#'
#' @examples
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

    num_newton <- as.integer(NA)
    num_grad <- result$num_grad
    if (is_implicit) {
      num_newton <- result$num_newton
    }

    create_onenode_tree(depth = 0,
                        invalid = result$error,
                        z = result$z1,
                        ham_system = ham_system,
                        num_newton = num_newton,
                        num_grad = num_grad,
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
  error <- NA
  if (H0-H1 >= 1000) {
    error <- "Divergence"
    warning("Divergence detected")
  }

  # return
  list(z1 = z1, error = error)
}

take_one_step_im <- function(z0, z_1, direction, ham_system) {

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

Norm <- function(v) {
  sqrt(sum(v^2))
}

newton_iteration <- function(g, Jg, x0, TOL = 1e-6, max_iter = 100) {

  eps <- TOL + 1
  x1 <- x0

  iter <- 0
  grads <- 0
  converged <- TRUE
  close_enough <- FALSE
  too_slow_counter <- 0
  # browser()
  # print("~~~~~~~~~~~~~~~~~~~~~~")
  eta <- 1e10
  prev_error <- c(0,0)
  prev_gmres_error <- c(0,0)
  quadratic_error <- c(0,0)
  while(eps > TOL) {
    g.x1 <- g(x1)
    # print(paste("Norm(g.new): ", Norm(g.x1)))
    # print(paste("grads: ", grads))
    # print("~~~~~~~~~~~~~~~~~~~~~~")
    # print(paste("x1:", x1))
    # print(paste("Norm(g.x1): ", Norm(g.x1)))
    # browser()



    # as.numeric(delta %*% -g.x1)
    # if (as.numeric(delta %*% -g.x1) < 0) {
    #   delta <- -g.x1
    # }

    #delta <- solve(diag(diag(Jg(x1))),-g.x1)
    # delta <- -g.x1

    # delta <- -g.x1
    # delta <- -as.numeric(t(Jg(x1)) %*% g.x1)
    # if (close_enough) {
    #   print("switched to newton")
    #   delta <- solve(Jg(x1),-g.x1)
    # }

    gmres_sol <- hessian_solve(x1, -g.x1, M, 2.0, rep(0.0, length(x1)), eta)
    delta <- gmres_sol$x
    prev_error <- Norm(g.x1)
    grads <- grads + gmres_sol$iterations

    # print("~~~~~~~~~~~~~~~~~~~~~~")
    # print(paste("-g.x1/Norm(-g.x1)", -g.x1/Norm(-g.x1)))
    # print(paste("delta/Norm(delta)", delta/Norm(delta)))
    # delta <- solve(Jg(x1),-g.x1)
    # jacobian_kept_delta <- all(sign(delta) == sign(-g.x1))
    # if (!jacobian_kept_delta) {
    #   delta <- -g.x1
    # }

    # print(paste("x1:", x1))
    # print(paste("delta:", delta))

    alpha <- 1
    line.search.iter <- 0
    g.new <- g.x1
    for(i in 1:20) {
      line.search.iter <- line.search.iter + 1

      # if delta is ridiculously big don't bother calculating it because that will
      # trigger a Stan error
      if(Norm(alpha*delta) > (1e1*Norm(x1))) {
        alpha <- alpha/2
        next
      }
      g.new <- g(x1 + alpha*delta)
      # print(paste("x1: ", x1))
      # print(paste("alpha*delta: ", alpha*delta))
      # print(paste("Norm(alpha*delta): ", Norm(alpha*delta)))
      # print(paste("(1e2*Norm(x1)): ", (1e2*Norm(x1))))
      # print(paste("g.new: ", g.new))
      if(Norm(g.new) >  Norm(g.x1)) {
        alpha <- alpha/2
      } else {
        break
      }
      #if(line.search.iter >= 9) browser()
    }
    # if (line.search.iter > 1) {
    #   warning("had to line search")
    # }
    # check if shrinking alpha actually reduces error further
    # for(i in 1:20) {
    #     line.search.iter <- line.search.iter + 1
    #
    #     alpha_new <- alpha/2.0
    #     g.new <- g(x1 + alpha*delta)
    #     g.new.new <- g(x1 + alpha_new*delta)
    #     if(Norm(g.new.new) <  Norm(g.new)) {
    #       alpha <- alpha_new
    #     } else {
    #       break
    #     }
    #     #if(line.search.iter >= 9) browser()
    # }

    if (Norm(g.new)/Norm(g.x1) >= 0.9) {
      too_slow_counter <- too_slow_counter + 1
    }
    if (too_slow_counter >= 3) {
      close_enough <- TRUE
    }
    if (iter >= 50) {
      close_enough <- TRUE
    }

    # print(paste("Norm(g.new): ", Norm(g.new)))
    # print(paste("Norm(g.new)/Norm(g.x1): ", Norm(g.new)/Norm(g.x1)))
    # print(paste("too_slow_counter: ", too_slow_counter))

    # this is for implementing 2.1 in eisenstat and walker
    prev_gmres_error <- g.x1 + as.vector(Jg(x1) %*% (alpha*delta))
    quadratic_error <- g(x1 + alpha*delta) - prev_gmres_error
    eta <- Norm(quadratic_error) / Norm(g.x1)

    grads <- grads + line.search.iter
    x1 <- x1 + alpha*delta
    #eps <- delta^2 %>% sum %>% sqrt
    eps <- Norm(g(x1))
    #print(eps)


    iter <- iter + 1
    if(iter >= max_iter) {
      warning("Newton Failed to Converge")
      converged <- FALSE
      break
    }
  }

  return(list(x1 = x1, num.newton.iters = iter, num_grad = grads, converged = converged))
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

leapfrog <- function(q0, p0, M, N, h) {

  D <- length(q0)

  q <- matrix(NA, nrow = N+1, ncol = D)
  p <- matrix(NA, nrow = N+1, ncol = D)
  omega <- rep(NA, N+1)

  q[1,] <- q0
  p[1,] <- p0
  omega[1] <- HessU(q[1,]) %>% eigen %>% .$values %>% max %>% sqrt

  # start solving
  for(n in 1:N) {

    phalf <- p[n,] - (h/2)*GradU(q[n,])
    q[n+1,] <- q[n,] + h*solve(M,phalf)
    p[n+1,] <- phalf - (h/2)*GradU(q[n+1,])
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
    mutate(H = H) %>%
    mutate(DeltaH = H0 - H) %>%
    mutate(acceptprob = pmin(1,exp(DeltaH)))

  return(soln)
}
