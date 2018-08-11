library(tidyverse)
library(mockery)
library(testthat)

source("StanHessianHelper/hessian_helper.R")
load_model("StanHessianHelper/neals_funnel.stan", list(dimension = 1))

M <- c(8.5, 56.3)
funnel <- create_hamiltonian_system(M, compute_U = U, compute_gradU = GradU)

q0 <- c(0,0)
p0 <- c(5.21,4.22)
h0 <- 0.5
lf <- create_integrator(is_implicit = FALSE)

mock_direction_sampler <- mock(c(1, -1, -1, -1,  1,  1,  1,  1, -1, -1), cycle = TRUE)
with_mock(sample = mock_direction_sampler, {
  sample <- get_single_nuts_sample(q0, p0, h0, ham_system = funnel, integrator = lf, max_treedepth = 10, DEBUG = TRUE)
})


