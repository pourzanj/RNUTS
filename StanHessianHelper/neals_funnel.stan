data {
  int<lower = 1> dimension;
}
parameters {
  real y;
  vector[dimension] x;
}
model {
  y ~ normal(0, 3);
  x ~ normal(0, exp(y/2));
}