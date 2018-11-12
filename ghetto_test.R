library(tidyverse)
library(mockery)
library(testthat)
library(rstan)

source("StanHessianHelper/hessian_helper.R")
load_model("StanHessianHelper/neals_funnel.stan", list(dimension = 1))

#########
load_model("StanHessianHelper/neals_funnel.stan", list(dimension = 100))

q0_1 <- rnorm(1,0, sd = 3)
q0 <- c(q0_1, rnorm(100, 0, sd = exp(q0_1/2)))
M <- c(8.5, rep(56.3, 100))
p0 <- sqrt(M) * rnorm(101)
impmid(q0, p0, diag(M), pos.first = TRUE, N = 1, h = 0.1)
####### END 100D
#############

M <- c(8.5, 56.3)
funnel <- create_hamiltonian_system(M, compute_U = U, compute_gradU = GradU)

q0 <- c(0,-1.5)
q0 <- c(-5,0)
p0 <- c(5.21,4.22)
lf <- create_integrator(is_implicit = FALSE)
im <- create_integrator(is_implicit = TRUE)

#### going down the funnel
true_soln <- leapfrog(c(0,0.01),c(-5,2),diag(M), N = 4096, h = 0.01) %>% select(-omega, -DeltaH) %>% gather(state,val,-t) %>% mutate(stepsize = "0.01")
lf_16 <- leapfrog(c(0,0.01),c(-5,2),diag(M), N = 256, h = 0.16) %>% select(-omega, -DeltaH) %>% gather(state,val,-t) %>% mutate(stepsize = "0.16")
lf_32 <- leapfrog(c(0,0.01),c(-5,2),diag(M), N = 128, h = 0.32) %>% select(-omega, -DeltaH) %>% gather(state,val,-t) %>% mutate(stepsize = "0.32")
im_32 <- impmid(c(0,0.01),c(-5,2),diag(M), N = 128, h = 0.32) %>% select(-omega, -DeltaH, -num.newton.iters, -num_grad, acceptprob = accept.prob) %>% gather(state,val,-t) %>% mutate(stepsize = "IM_0.32")
im_64 <- impmid(c(0,0.01),c(-5,2),diag(M), N = 64, h = 0.64) %>% select(-omega, -DeltaH, -num.newton.iters, -num_grad, acceptprob = accept.prob) %>% gather(state,val,-t) %>% mutate(stepsize = "IM_0.64")
im_128 <- impmid(c(0,0.01),c(-5,2),diag(M), N = 32, h = 1.28) %>% select(-omega, -DeltaH, -num.newton.iters, -num_grad, acceptprob = accept.prob) %>% gather(state,val,-t) %>% mutate(stepsize = "IM_1.28")
im_256 <- impmid(c(0,0.01),c(-5,2),diag(M), N = 32, h = 2.56) %>% select(-omega, -DeltaH, -num.newton.iters, -num_grad, acceptprob = accept.prob) %>% gather(state,val,-t) %>% mutate(stepsize = "IM_2.56")


numerical_solns <- bind_rows(im_64, im_128) %>% mutate(stepsize = as.factor(stepsize))
true_soln %>% ggplot(aes(t, val)) + geom_line() + facet_grid(state ~ ., scales = "free") +
  geom_line(aes(color = factor(stepsize)), data = numerical_solns) + geom_point(aes(color = factor(stepsize)),data = numerical_solns)

true_soln %>% ggplot(aes(t, val)) + geom_line() + facet_grid(state ~ ., scales = "free")
######

q0 <- funnel.draws %>% sample_n(1) %>% select(q1,q2) %>% as.matrix %>% as.vector()
p0 <- funnel$get_momentum_sample() %>% as.vector
sample <- get_single_nuts_sample(q0, p0, h0 = 1.0, ham_system = funnel, integrator = im, max_treedepth = 10, DEBUG = TRUE)
sample$hist

root_sample <- sample$hist %>% filter(is.na(depth))
root_sample_q0 <- root_sample %>% select(q1,q2) %>% as.matrix %>% as.vector
root_sample_p0 <- root_sample %>% select(p1,p2) %>% as.matrix %>% as.vector
NA_index <- which(is.na(sample$hist$depth))

backwards_traj <- impmid(root_sample_q0, -root_sample_p0, diag(M), pos.first = TRUE, N = (NA_index-1+60), h = 1.0) %>% mutate(t = row_number()) %>% arrange(desc(t)) %>% filter(t > 1) %>% mutate(p1 = -p1, p2 = -p2)
forwards_traj <- impmid(root_sample_q0, root_sample_p0, diag(M), pos.first = TRUE, N = (nrow(sample$hist)-NA_index+60), h = 1.0)

bind_rows(backwards_traj, forwards_traj) %>%
  mutate(t = 1:nrow(.)) %>%
  select(-DeltaH,-num.newton.iters,-omega, -num_grad) %>%
  gather(state,val,-t) %>%
  ggplot(aes(t,val)) +
  geom_point() +
  geom_line() +
  facet_grid(state ~ ., scales = "free") +
  geom_point(color = "red", data = sample$hist %>% select(-depth,-invalid,-num_newton,-num_grad) %>% mutate(t = 61:(61+nrow(sample$hist)-1)) %>% gather(state,val,-t)) +
  geom_line(color = "red", data = sample$hist %>% select(-depth,-invalid,-num_newton,-num_grad) %>% mutate(t = 61:(61+nrow(sample$hist)-1)) %>% gather(state,val,-t))

PlotFunnelNutsDraw(GenerateFunnelPlot(-2,2,min(sample$hist$q1)-1,max(sample$hist$q1)+1), sample$hist, sample$q$q)
PlotFunnelNutsDraw(GenerateFunnelPlot(-2,2,min(sample$hist$q1)-1,max(sample$hist$q1)+1), sample$hist, sample$q$q)

sample$hist %>%
  select(-invalid,-num_newton,-num_grad) %>%
  mutate(t = seq(0,2*nrow(.)-1,by=2)) %>%
  gather(state,val,-t,-depth) %>%
  mutate(depth=factor(depth)) %>%
  ggplot(aes(t,val)) +
  geom_line() +
  geom_point(aes(color = depth)) +
  facet_grid(state ~ ., scales = "free")

mock_direction_sampler <- mock(c(-1, -1, -1, -1,  1,  1,  1,  1, -1, -1), cycle = TRUE)
with_mock(sample = mock_direction_sampler, {
  sample <- get_single_nuts_sample(q0, p0, h0 = 2.0, ham_system = funnel, integrator = lf, max_treedepth = 10, DEBUG = TRUE)
})

mock_direction_sampler <- mock(c(-1, -1, -1, -1,  1,  1,  1,  1, -1, -1), cycle = TRUE)
with_mock(sample = mock_direction_sampler, {
  sample <- get_single_nuts_sample(q0, p0, h0 = 2.0, ham_system = funnel, integrator = im, max_treedepth = 10, DEBUG = TRUE)
})
sample$hist
PlotFunnelNutsDraw(GenerateFunnelPlot(-2,2,min(sample$hist$q1)-1,max(sample$hist$q1)+1), sample$hist, sample$q$q)

N <- 1000
funnel.draws <- tibble(q1 = rnorm(N, 0, 3)) %>%
  mutate(q2 = rnorm(1000,0,exp(q1/2)))


Uvec <- Vectorize(function(q1, q2) U(c(q1,q2)))
GenerateFunnelPlot <- function(xmin,xmax,ymin,ymax) {

  expand.grid(q2 = seq(xmin,xmax,(xmax-xmin)/200), q1 = seq(ymin,ymax,(ymax-ymin)/200)) %>%
    as_tibble %>%
    mutate(U = Uvec(q1,q2)) %>%
    ggplot(aes(q2,q1)) +
    geom_raster(aes(fill = log(U+2))) +
    scale_fill_gradientn(colours = c("white", "orange", "black"))

}

Uvec <- Vectorize(function(q1, q2) U(c(q1,q2)))
GenerateFunnelSurface <- function(xmin,xmax,ymin,ymax) {

  q2_grid <- seq(xmin,xmax,(xmax-xmin)/200)
  q1_grid <- seq(ymin,ymax,(ymax-ymin)/200)
  U <- outer(q1_grid, q2_grid, FUN = Gvec)

  plot_ly(x = q2_grid, y = q1_grid, z = U) %>% add_surface()

}

p + geom_point(color = "red", data = tibble(q1 = x0[1], q2 = x0[2])) +
  geom_segment(aes(xend = q2 + eps*delta_grad[2], yend = q1 + eps*delta_grad[1]),arrow = arrow(length = unit(0.1,"cm")), color = "red", data = tibble(q1 = x0[1], q2 = x0[2])) +
  geom_segment(aes(xend = q2 + eps2*delta_newton[2], yend = q1 + eps2*delta_newton[1]),arrow = arrow(length = unit(0.1,"cm")), color = "green", data = tibble(q1 = x0[1], q2 = x0[2])) +
  geom_segment(aes(xend = q2 + eps3*delta_vanilla[2], yend = q1 + eps3*delta_vanilla[1]),arrow = arrow(length = unit(0.1,"cm")), color = "blue", data = tibble(q1 = x0[1], q2 = x0[2])) +
  geom_segment(aes(xend = q2 + eps4*delta_jeng[2], yend = q1 + eps4*delta_jeng[1]),arrow = arrow(length = unit(0.1,"cm")), color = "yellow", data = tibble(q1 = x0[1], q2 = x0[2])) +
  geom_segment(aes(xend = q2 + eps5*delta_abs[2], yend = q1 + eps5*delta_abs[1]),arrow = arrow(length = unit(0.1,"cm")), color = "orange", data = tibble(q1 = x0[1], q2 = x0[2]))

funnel.plot <- GenerateFunnelPlot(-70,70,-10,10)
funnel.plot
funnel.plot + geom_point(data = funnel.draws) + xlim(-70,70) + ylim(-10,10)

PlotFunnelNutsDraw <- function(funnel.plot, nuts.draw.hist, q) {

  origin <- nuts.draw.hist %>% filter(is.na(depth))
  left.most.point <- head(nuts.draw.hist,1)
  right.most.point <- tail(nuts.draw.hist,1)

  left.vector <- tibble(q0.1 = origin$q1, q0.2 = origin$q2, q1.1 = left.most.point$q1, q1.2 = left.most.point$q2)
  right.vector <- tibble(q0.1 = origin$q1, q0.2 = origin$q2, q1.1 = right.most.point$q1, q1.2 = right.most.point$q2)
  uturn <- nuts.draw.hist %>% filter(invalid == "U-Turn")
  head.uturn <- head(uturn,1)
  tail.uturn <- tail(uturn,1)

  funnel.plot +
    #geom_point(data = funnel.draws %>% filter(q2 < 100 & q2 > -100)) +
    geom_path(data = nuts.draw.hist) +
    geom_point(size = 1.7, data = nuts.draw.hist %>% mutate(depth = as.integer(depth))) +
    geom_point(aes(color=as.factor(depth)), size = 1.2,
               data = nuts.draw.hist %>% mutate(depth = as.integer(depth))) +
    geom_point(color = "red", size = 3.0, data = nuts.draw.hist %>% filter(is.na(depth))) +
    geom_point(color = "green",  size = 3.0,data = tibble(q1 = q[1], q2 = q[2]))
  # +
    # geom_segment(aes(xend = q2+p2, yend = q1+p1), color = "red",
    #              arrow = arrow(length = unit(0.01, "npc")),
    #              nuts.draw.hist %>% filter(is.na(depth)))
}

PlotFunnelNutsDraw(GenerateFunnelPlot(-0.5,0.5,-10,-2.5), sample$hist)

impmid(q0, p0, diag(M), pos.first = TRUE, N = 10, h = 2) %>%
  select(-DeltaH,-num.newton.iters,-omega) %>%
  gather(state,val,-t) %>%
  ggplot(aes(t,val)) +
  geom_point() +
  geom_line() +
  facet_grid(state ~ ., scales = "free")

system.time(samples_im <- get_nuts_samples(num_samples = 1000, q0 = q0, h0 = 2.0, ham_system = funnel, integrator = im, DEBUG = FALSE))
system.time(samples_lf <- get_nuts_samples(num_samples = 5, q0 = q0, h0 = 0.01, ham_system = funnel, integrator = lf, max_treedepth = 12,  DEBUG = TRUE))


effectiveSize(samples_im$q1)
effectiveSize(samples_lf$q1)

effectiveSize(samples_im$q2)
effectiveSize(samples_lf$q2)
