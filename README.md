# RNUTS

This repo contains all the code for our planned paper on Stiffness and using Implicit Midpoint in NUTS.

# High-Level Goals of Paper
- teach stats community about stiffness and how explicit integrator fails (hit by section 1 of paper)
- show with examples this is important in practice because it can speed up sampling (hit by 3)
- show which example explicit integrator is more prone to fail on (hit by 1 a little when we define stiffness, hit by 3)
- get Stan community interested enough to want to have midpoint and new warmup implemented in production (4)
- eventually need to make decusion of whether paper should just be about midpoint or include warmup

## R Code to do
- get Ben's GMRES incorporated in the code
- clean up code? ask Dan for advice on this. is it worth it?
- should we do unit tests?
- more documentation

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
- 250-D Gaussian. Do a plot of Work per effective sample on y-axis vs. largest correlation coefficent on x-axis. Leapfrog will be one line and Midpoint will be one line. (maybe a 2D gaussian would be preferred)
- possibly point out an example where non-clustered eigenvalues can take longer
- possible example of linear regression w/o QR or with non-centered covariates which leads to high posterior correlation
- Neal's Funnel. popular example and illustrates how one has clustered eigenvalues when doing a hierarchical model
- Latent variable e.g. a GP with non-gaussian noise. the points of the GP will have high posterior correlation (Bales)
- PKPD with high posterior correlation (gelman said to include a PKPD since it's hot)
- estimating a covariance matrix using the Cholesky factor parameterization is highly stiff but non-Gaussian

4. New proposed warmup for HMC
    - check stiff periodically which is really cheap using power iterations
    - determine whether to use Implicit or Explicit solver
    - if we use Explicit solver use the stepsize chosen by hessian eigenvalue NOT from nesterov optimization. this way we'll be less likely to have divergences

## Things to Discuss
- getting in to Stan (Dan, Bob, and Bales)
- how to test for correctness? just look at quantiles of posterior converging?
- should numerical experiments showing correctness of samples be included in the paper? possibly in supplemtary material. showing that it's symplectic is enough
- how should examples be done? my R code leapfrog vs midpoint
- use stepsize and mass matrix Stan decides? TBD how to determine midpoint stepsize. show a plot of work/sample vs. stepsize
