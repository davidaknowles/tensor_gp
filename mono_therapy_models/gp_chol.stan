data {
  int<lower=0> N;
  matrix[N,N] sqDist; 
  vector[N] y;
}
parameters {
  real<lower=0> eta_sq;
  real<lower=0> inv_rho_sq;
  real<lower=0> sigma_sq;
  vector[N] z;
}
transformed parameters {
      real<lower=0> rho_sq;
      rho_sq <- inv(inv_rho_sq);
}
model {
  matrix[N,N] Sigma;
  matrix[N,N] L;
  // off-diagonal elements
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      Sigma[i,j] <- eta_sq * exp(-rho_sq * sqDist[i,j]);
      Sigma[j,i] <- Sigma[i,j];
    } 
  }
  for (k in 1:N)
    Sigma[k,k] <- eta_sq + sigma_sq;
  
  L <- cholesky_decompose(Sigma);
  
  eta_sq ~ cauchy(0,5);
  inv_rho_sq ~ cauchy(0,5);
  sigma_sq ~ cauchy(0,5);
  z ~ normal(0,1);
  y <- L * z; 
}