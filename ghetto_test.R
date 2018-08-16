library(tidyverse)
library(mockery)
library(testthat)

source("StanHessianHelper/hessian_helper.R")
load_model("StanHessianHelper/neals_funnel.stan", list(dimension = 1))

M <- c(8.5, 56.3)
funnel <- create_hamiltonian_system(M, compute_U = U, compute_gradU = GradU)

q0 <- c(0,-1.5)
q0 <- c(-5,0)
p0 <- c(5.21,4.22)
lf <- create_integrator(is_implicit = FALSE)
im <- create_integrator(is_implicit = TRUE)

q0 <- funnel.draws %>% sample_n(1) %>% select(q1,q2) %>% as.matrix %>% as.vector()
p0 <- funnel$get_momentum_sample() %>% as.vector
sample <- get_single_nuts_sample(q0, p0, h0 = 2.0, ham_system = funnel, integrator = lf, max_treedepth = 10, DEBUG = TRUE)
sample$hist

head_sample <- head(sample$hist,1)
impmid(head_sample %>% select(q1,q2) %>% as.matrix %>% as.vector,
       head_sample %>% select(p1,p2) %>% as.matrix %>% as.vector,
       diag(M), pos.first = TRUE, N = 13, h = 2.0) %>%
  impmid(q00,p00,diag(M),TRUE,N=8,h=2.0) %>%
  select(-DeltaH,-num.newton.iters,-omega) %>%
  gather(state,val,-t) %>%
  ggplot(aes(t,val)) +
  geom_point() +
  geom_line() +
  facet_grid(state ~ ., scales = "free") +
  geom_point(color = "red", data = sample$hist %>% select(-depth,-invalid,-num_newton,-num_grad) %>% mutate(t = seq(0,2*nrow(.)-1,by=2)) %>% gather(state,val,-t)) +
  geom_line(color = "red", data = sample$hist %>% select(-depth,-invalid,-num_newton,-num_grad) %>% mutate(t = seq(0,2*nrow(.)-1,by=2)) %>% gather(state,val,-t))

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

p + geom_point(color = "red", data = tibble(q1 = x0[1], q2 = x0[2])) + geom_segment(aes(xend = q2 + eps*delta_grad[2], yend = q1 + eps*delta_grad[1]),arrow = arrow(length = unit(0.1,"cm")), color = "red", data = tibble(q1 = x0[1], q2 = x0[2])) + geom_segment(aes(xend = q2 + eps*delta_newton[2], yend = q1 + eps*delta_newton[1]),arrow = arrow(length = unit(0.1,"cm")), color = "red", data = tibble(q1 = x0[1], q2 = x0[2]))

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

samples_im <- get_nuts_samples(num_samples = 1000, q0 = q0, h0 = 2.0, ham_system = funnel, integrator = im, DEBUG = TRUE)
samples_lf <- get_nuts_samples(num_samples = 1000, q0 = q0, h0 = 2.0, ham_system = funnel, integrator = lf, DEBUG = TRUE)
