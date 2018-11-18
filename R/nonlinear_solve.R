#' Create generalized Newton solve function.
#'
#' A function for creating a custom function that solves g(x)=0. Special cases include Newton iterations and Newton-Krylov.
#'
#' This function returns a function to solve nonlinear equations using a given select_direction,
#' and select_stepsize function. Internally, the function will be updating the intial guess as
#' x = x + alpha*delta.
#'
#' The select_direction function selects the vector delta. For Newton iterations
#' this would simply be delta = solve(Jg(x), -g(x)).
#'
#' The select_stepsize function selects alpha. For plain Newton iterations this would always be
#' one, but it can be replaced by more elaborate techniques such as line search.
#'
#' @param select_direction
#' @param select_stepsize
#' @param TOL
#' @param MAX_ITER
#'
#' @return
#' @export
#'
#' @examples
create_generalized_newton_solve <- function(select_direction, select_stepsize, compute_residual, TOL = 1e-6, MAX_ITER = 100) {

  do_generalized_newton_solve <- function(g, Jg, Jg_v, x0) {

    eps <- TOL + 1
    x <- x0

    iter <- 0
    num_hess_evals <- 0
    num_hess_vec_prod_evals <- 0
    converged <- TRUE
    prev_error <- c(0,0)

    while(eps > TOL) {
      g_x <- g(x)

      # 1) choose direction to update, delta
      direction_soln <- select_direction(g_x, Jg, Jg_v)
      delta <- direction_soln$x
      num_hess_evals <- num_hess_evals + direction_soln$num_hess_evals
      num_hess_vec_prod_evals <- num_hess_vec_prod_evals + direction_soln$num_hess_vec_prod_evals

      # 2) choose stepsize
      stepsize_soln <- select_direction(g_x, Jg, Jg_v)
      alpha <- stepsize_soln$x
      num_hess_evals <- num_hess_evals + stepsize_soln$num_hess_evals
      num_hess_vec_prod_evals <- num_hess_vec_prod_evals + stepsize_soln$num_hess_vec_prod_evals

      # 3) make update
      x <- x + alpha*delta

      # 4) check convergence

      eps <- compute_residual()

      iter <- iter + 1
      if(converged) {

        if(iter >= MAX_ITER) {
          warning("Newton reached max number of iterations.")
          did_converge <- FALSE
          break
        }
      }


    }

    return(list(x1 = x1, did_converge = converged, num_newton_iters = iter,
                num_hess_evals = grads, num_hess_vec_prod_evals))
  }

}



# generalized_newton_solve <- function(g, Jg, x0, TOL = 1e-6, max_iter = 100) {
#
#   eps <- TOL + 1
#   x1 <- x0
#
#   iter <- 0
#   grads <- 0
#   converged <- TRUE
#   close_enough <- FALSE
#   too_slow_counter <- 0
#   # browser()
#   # print("~~~~~~~~~~~~~~~~~~~~~~")
#   eta <- 1e10
#   prev_error <- c(0,0)
#   prev_gmres_error <- c(0,0)
#   quadratic_error <- c(0,0)
#   while(eps > TOL) {
#     g.x1 <- g(x1)
#     # print(paste("Norm(g.new): ", Norm(g.x1)))
#     # print(paste("grads: ", grads))
#     # print("~~~~~~~~~~~~~~~~~~~~~~")
#     # print(paste("x1:", x1))
#     # print(paste("Norm(g.x1): ", Norm(g.x1)))
#     # browser()
#
#
#
#     # as.numeric(delta %*% -g.x1)
#     # if (as.numeric(delta %*% -g.x1) < 0) {
#     #   delta <- -g.x1
#     # }
#
#     #delta <- solve(diag(diag(Jg(x1))),-g.x1)
#     # delta <- -g.x1
#
#     # delta <- -g.x1
#     # delta <- -as.numeric(t(Jg(x1)) %*% g.x1)
#     # if (close_enough) {
#     #   print("switched to newton")
#     #   delta <- solve(Jg(x1),-g.x1)
#     # }
#
#     gmres_sol <- hessian_solve(x1, -g.x1, M, 2.0, rep(0.0, length(x1)), eta)
#     delta <- gmres_sol$x
#     prev_error <- Norm(g.x1)
#     grads <- grads + gmres_sol$iterations
#
#     # print("~~~~~~~~~~~~~~~~~~~~~~")
#     # print(paste("-g.x1/Norm(-g.x1)", -g.x1/Norm(-g.x1)))
#     # print(paste("delta/Norm(delta)", delta/Norm(delta)))
#     # delta <- solve(Jg(x1),-g.x1)
#     # jacobian_kept_delta <- all(sign(delta) == sign(-g.x1))
#     # if (!jacobian_kept_delta) {
#     #   delta <- -g.x1
#     # }
#
#     # print(paste("x1:", x1))
#     # print(paste("delta:", delta))
#
#     alpha <- 1
#     line.search.iter <- 0
#     g.new <- g.x1
#     for(i in 1:20) {
#       line.search.iter <- line.search.iter + 1
#
#       # if delta is ridiculously big don't bother calculating it because that will
#       # trigger a Stan error
#       if(Norm(alpha*delta) > (1e1*Norm(x1))) {
#         alpha <- alpha/2
#         next
#       }
#       g.new <- g(x1 + alpha*delta)
#       # print(paste("x1: ", x1))
#       # print(paste("alpha*delta: ", alpha*delta))
#       # print(paste("Norm(alpha*delta): ", Norm(alpha*delta)))
#       # print(paste("(1e2*Norm(x1)): ", (1e2*Norm(x1))))
#       # print(paste("g.new: ", g.new))
#       if(Norm(g.new) >  Norm(g.x1)) {
#         alpha <- alpha/2
#       } else {
#         break
#       }
#       #if(line.search.iter >= 9) browser()
#     }
#     # if (line.search.iter > 1) {
#     #   warning("had to line search")
#     # }
#     # check if shrinking alpha actually reduces error further
#     # for(i in 1:20) {
#     #     line.search.iter <- line.search.iter + 1
#     #
#     #     alpha_new <- alpha/2.0
#     #     g.new <- g(x1 + alpha*delta)
#     #     g.new.new <- g(x1 + alpha_new*delta)
#     #     if(Norm(g.new.new) <  Norm(g.new)) {
#     #       alpha <- alpha_new
#     #     } else {
#     #       break
#     #     }
#     #     #if(line.search.iter >= 9) browser()
#     # }
#
#     if (Norm(g.new)/Norm(g.x1) >= 0.9) {
#       too_slow_counter <- too_slow_counter + 1
#     }
#     if (too_slow_counter >= 3) {
#       close_enough <- TRUE
#     }
#     if (iter >= 50) {
#       close_enough <- TRUE
#     }
#
#     # print(paste("Norm(g.new): ", Norm(g.new)))
#     # print(paste("Norm(g.new)/Norm(g.x1): ", Norm(g.new)/Norm(g.x1)))
#     # print(paste("too_slow_counter: ", too_slow_counter))
#
#     # this is for implementing 2.1 in eisenstat and walker
#     prev_gmres_error <- g.x1 + as.vector(Jg(x1) %*% (alpha*delta))
#     quadratic_error <- g(x1 + alpha*delta) - prev_gmres_error
#     eta <- Norm(quadratic_error) / Norm(g.x1)
#
#     grads <- grads + line.search.iter
#     x1 <- x1 + alpha*delta
#     #eps <- delta^2 %>% sum %>% sqrt
#     eps <- Norm(g(x1))
#     #print(eps)
#
#
#     iter <- iter + 1
#     if(iter >= max_iter) {
#       warning("Newton Failed to Converge")
#       converged <- FALSE
#       break
#     }
#   }
#
#   return(list(x1 = x1, num.newton.iters = iter, num_grad = grads, converged = converged))
# }
