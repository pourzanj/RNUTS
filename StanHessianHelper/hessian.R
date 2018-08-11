library(tidyverse)
library(ggplot2)
library(rstan)
library(Rcpp)

load_model = function(model_file, params_list) {
  envir = list2env(params_list)
  if(length(params_list) == 0) {
    system("touch data.dump")
  } else {
    stan_rdump(names(params_list), file = "data.dump", envir = envir)
  }
  
  if(system(paste0('/Applications/cmdstan-2.17.1/bin/stanc --allow_undefined --name=linear_regression_model ', model_file)) != 0) {
    stop("model '", model_file,"' does not exist")
  }
  system('mv linear_regression_model.cpp StanHessianHelper/linear_regression_model.hpp')
  system(paste0('cp ', paste0(getSrcDirectory(function(dummy) {dummy})), 'StanHessianHelper/helper.cpp helper.cpp'))
  sourceCpp("StanHessianHelper/helper.cpp")
  set_data("data.dump")
}

a = 1.5
b = 0.5
sigma = 0.25
N = 20

x = seq(0.0, 1.0, length = N)
(y = sapply(x, function(xv) rnorm(1, a * xv + b, sigma)))

list(x = x, y = y) %>% as.tibble %>%
  ggplot(aes(x, y)) +
  geom_point()

load_model("StanHessianHelper/linear_regression.stan", list(N = N, x = x, y = y))

jacobian(c(a, b, sigma))
h = hessian(c(a, b, sigma))
vec = rnorm(3)
hv = hessian_vector(c(a, b, sigma), vec)
dt = 2.0
M = abs(rnorm(3))
solve(diag(rep(1, 3)) + (dt * dt / 4.0) * diag(1 / M) %*% h$hess, vec)
hessian_solve(c(a, b, sigma), vec, M, dt, rep(0, 3), 1e-10)$x
eigen(diag(rep(1, 3)) + (dt * dt / 4.0) * diag(1 / M) %*% h$hess)

model = stan_model("models/linear_regression.stan")
fit = optimizing(model, data = list(N = N, x = x, y = y), hessian = TRUE)

(jac = jacobian(fit$par[1:3]))
(hess = hessian(fit$par[1:3]))

samples = sampling(model, data = list(N = N, x = x, y = y), chains = 1)

pairs(samples, pars = c("a", "b", "sigma_log"))

eigen(hess)