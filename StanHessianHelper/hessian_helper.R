library(rstan)
library(Rcpp)

load_model <- function(model_file, params_list) {
  envir = list2env(params_list)
  if(length(params_list) == 0) {
    system("touch data.dump")
  } else {
    stan_rdump(names(params_list), file = "data.dump", envir = envir)
  }
  
  # create linear_regression_model.cpp from Stan file
  if(system(paste0('/Applications/cmdstan-2.17.1/bin/stanc --allow_undefined --name=linear_regression_model ', model_file)) != 0) {
    stop("model '", model_file,"' does not exist")
  }
  
  # move file to StanHessianHelper and rename it with .hpp
  #system(paste0('mv linear_regression_model.cpp ', system.file(package = "RNUTS"), '/StanHessianHelper/linear_regression_model.hpp'))
  system(paste0('mv linear_regression_model.cpp ', 'StanHessianHelper/linear_regression_model.hpp'))
  
  # compile with RCPP
  # system(paste0('cp ', paste0(getSrcDirectory(function(dummy) {dummy})), 'helper.cpp helper.cpp'))
  #sourceCpp(paste0(system.file(package = "RNUTS"), "/StanHessianHelper/helper.cpp"))
  sourceCpp(paste0("StanHessianHelper/helper.cpp"))
  set_data("data.dump")

  U <<- function(q) -jacobian(q)$u
  GradU <<- function(q) -jacobian(q)$jac
  HessU <<- function(q) -hessian(q)$hess
  HessVecProd <<- function(q, vec) -hessian_vector(q, vec)$hessv
}

#load_model("nuts_tests/neals_funnel/neals_funnel.stan", list(dimension = 1))
