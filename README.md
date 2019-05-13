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

### Adding Tests
- since code is still kind of unstable it may only be worth it to run tests for outside code such as the GMRES solver
- have a top level script that runs all tests and make sure it always works

## Things to Discuss
- getting in to Stan
- how to test for correctness? just look at quantiles of posterior converging?
- should numerical experiments showing correctness of samples be included in the paper? possibly in supplemtary material. showing that it's symplectic is enough
- use stepsize and mass matrix Stan decides? TBD how to determine midpoint stepsize. show a plot of work/sample vs. stepsize
