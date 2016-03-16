data {
  int<lower=0> N;
  int<lower=0> Ntest;
  matrix[N,N] sqDist; 
  vector[N] y;
  matrix[Ntest,Ntest] sqDistTest; 
  matrix[N,Ntest] sqDistTrainTest; 
}
transformed data {
  vector[N] mu;
  for (i in 1:N) mu[i] <- 0;
}
parameters {
  real<lower=0> eta_sq;
  real<lower=0> inv_rho_sq;
  real<lower=0> sigma_sq;
}
transformed parameters {
  matrix[N,N] Sigma;
  real<lower=0> rho_sq;
  rho_sq <- inv(inv_rho_sq);
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      Sigma[i,j] <- eta_sq * exp(-rho_sq * sqDist[i,j]);
      Sigma[j,i] <- Sigma[i,j];
    } 
  }
  for (k in 1:N)
    Sigma[k,k] <- eta_sq + sigma_sq;
}
model {

  eta_sq ~ cauchy(0,5);
  inv_rho_sq ~ cauchy(0,5);
  sigma_sq ~ cauchy(0,5);
  y ~ multi_normal(mu,Sigma);
}
generated quantities {
  vector[Ntest] testMean;
  matrix[Ntest,Ntest] Tau;
  {
    matrix[Ntest,Ntest] Omega;
    matrix[N,Ntest] K;
    matrix[Ntest,N] K_transpose_div_Sigma;
    for (i in 1:Ntest)
      for (j in 1:Ntest)
        Omega[i,j] <- eta_sq * exp(-rho_sq * sqDistTest[i,j]) + if_else(i==j, sigma_sq, 0.0);
    for (i in 1:N)
      for (j in 1:Ntest)
        K[i,j] <- eta_sq * exp(-rho_sq * sqDistTrainTest[i,j]); 
    K_transpose_div_Sigma <- K' / Sigma;
    testMean <- K_transpose_div_Sigma * y;
    Tau <- Omega - K_transpose_div_Sigma * K;
    for (i in 1:(Ntest-1))
      for (j in (i+1):Ntest)
        Tau[i,j] <- Tau[j,i];
  }
}