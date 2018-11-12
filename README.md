# RNUTS

This repo contains all the code for our planned paper on Stiffness and using Implicit Midpoint in NUTS.

# Outline of Code
get_nuts_samples(z0) # nuts.R
   1) ham_system <- create_hamiltonian_system() # nuts.R
   2) integrator <- create_integrator() # integrators.R
   3) get_single_nuts_sample(z0, ham_system, integrator) #nuts.R
      A) p0 <- ham_system$get_momentum_sample()
      B) tree <- create_onenode_tree(z0, ham_system) # initialize tree from build_tree.R
      C) new_subtree <- build_tree(depth, tree$z_plus, tree$z_plus_1, directions[depth+1], ham_system, integrator, DEBUG)

# High-Level Goals of Paper
- teach stats community about stiffness and how explicit integrator fails (hit by section 1 of paper)
- show with examples this is important in practice because it can speed up sampling (hit by 3)
- show which example explicit integrator is more prone to fail on (hit by 1 a little when we define stiffness, hit by 3)
- get Stan community interested enough to want to have midpoint and new warmup implemented in production (4)
- eventually need to make decusion of whether paper should just be about midpoint or include warmup

### Adding Tests
- since code is still kind of unstable it may only be worth it to run tests for outside code such as the GMRES solver
- have a top level script that runs all tests and make sure it always works

## Rough Overview of Paper
1. Talk about stability and stiffness.
- How highly correlated posteriors lead to stiff problems which can be measured by largest eigenvalue of hessian
- What is stiffness and how explicit integrators go unstable (divergence)
- if you have a stiff problem an explicit integrator will be forced to take tiny steps which is computationally inefficient
- explain Implicit Midpoint and why doesn't go unstable
- since it's symplectic we can sub it in for leapfrog in HMC (explain a little what symplectic is and why it's valid)

2. Practical use of Implicit Midpoint. NEWTON IS OK BUT NEWTON-KRYLOV IS THE WAY TO GO
- to do an implicit integrator we need a solve
- to do a solve efficiently we can't just use fixed point iteration we have to use newton
- to use newton a hessian and hessians scale quadratically with the dimension of parameter space. 
- isntead of vanilla Newton iterations we can use Newton Krylov which only require hessian-vector products which scale linearly in problem size
- point out how newton krylov worst case can happen eigenvalues of hessian are spread out, but in bayesian models eigenvalues are usually clustered because parameters are in some hierarchy
- other bells and whistles, line search, GMRES forcing term (mention in passing)

3. Examples in HMC
GAUSSIAN EXAMPLES
a) 2D Gaussian. Do a plot of Work per effective sample on y-axis vs. largest correlation coefficent on x-axis. Leapfrog will be one line and Midpoint will be one line
b) Aki temperature linear regression w/o non-centered covariates which leads to high posterior correlation
c) 250-D Guassian point out how non-clustered eigenvalues can take longer because of Newton-Krylov convergence properties
NON-GAUSSIAN EXAMPLEs
d) Neal's Funnel. popular example and illustrates how one has clustered eigenvalues when doing a hierarchical model
e) Latent variable e.g. a GP with non-gaussian noise. the points of the GP will have high posterior correlation (Bales)
f) PKPD with high posterior correlation

4. New proposed warmup for HMC
    - check stiff periodically which is really cheap using power iterations
    - determine whether to use Implicit or Explicit solver
    - if we use Explicit solver use the stepsize chosen by hessian eigenvalue NOT from nesterov optimization. this way we'll be less likely to have divergences

## Things to Discuss
- getting in to Stan
- how to test for correctness? just look at quantiles of posterior converging?
- should numerical experiments showing correctness of samples be included in the paper? possibly in supplemtary material. showing that it's symplectic is enough
- use stepsize and mass matrix Stan decides? TBD how to determine midpoint stepsize. show a plot of work/sample vs. stepsize
