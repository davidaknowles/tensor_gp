data {
  int<lower=0> K; 
  int<lower=0> D; 
  int<lower=0> C;
  int<lower=0> N;
  int<lower=0,upper=C> drug[N]; 
  int<lower=0,upper=C> cellLine[N]; 
  int<lower=0> Ntest;
  int<lower=0,upper=C> drugtest[Ntest]; 
  int<lower=0,upper=C> cellLinetest[Ntest];
  matrix[C,C] xx; 
  real synergy[N];
}
transformed data {
  vector[C] mu;
  for (i in 1:C) mu[i] <- 0;
}
parameters {
  real drugV[D,K];
  vector[C] z[K];
  real<lower=0> noise; 
  real meanSyn;
  real<lower=0> drugscale[K];
  real<lower=0> scalescale;
  real<lower=0> eta_sq;
  real<lower=0> sigma_sq;
}
transformed parameters {
  vector[C] cellV[K];
  matrix[C,C] Sigma;
  matrix[C,C] L;
  for (i in 1:(C-1)) {
    for (j in (i+1):C) {
      Sigma[i,j] <- eta_sq * xx[i,j];
      Sigma[j,i] <- Sigma[i,j];
    } 
  }
  for (k in 1:C)
    Sigma[k,k] <- eta_sq + sigma_sq;
  L <- cholesky_decompose(Sigma);
  for (k in 1:K) 
    cellV[k] <- L * z[k]; 
}
model {
  eta_sq ~ lognormal(0,1);
  sigma_sq ~ lognormal(0,1);

  scalescale ~ cauchy(0,1);
  drugscale ~ cauchy(0,scalescale);
  for (d in 1:D)
    for (k in 1:K)
      drugV[d,k] ~ normal(0,drugscale[k]);
  
  noise ~ cauchy(0,5);

  for (n in 1:N) {
    real temp[K]; 
    for (k in 1:K)
      temp[k] <- drugV[drug[n],k] * cellV[k][cellLine[n] ]; 
    synergy[n] ~ normal( meanSyn + sum(temp), noise );
  }
}
generated quantities {
  real test_synergy[Ntest];
  for (n in 1:Ntest) {
    real temp[K]; 
    for (k in 1:K)
      temp[k] <- drugV[drugtest[n],k] * cellV[k][cellLinetest[n] ]; 
    test_synergy[n] <- meanSyn + sum(temp);
  }
}