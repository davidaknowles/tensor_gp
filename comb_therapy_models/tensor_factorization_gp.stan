data {
  int<lower=0> K; 
  int<lower=0> D; 
  int<lower=0> C;
  int<lower=0> N;
  int<lower=0,upper=C> drugA[N]; 
  int<lower=0,upper=C> drugB[N];
  int<lower=0,upper=C> cellLine[N]; 
  int<lower=0> Ntest;
  int<lower=0,upper=C> drugAtest[Ntest]; 
  int<lower=0,upper=C> drugBtest[Ntest];
  int<lower=0,upper=C> cellLinetest[Ntest];
  matrix[C,C] sqDist; 
  real synergy[N];
}
transformed data {
  vector[C] mu;
  for (i in 1:C) mu[i] <- 0;
}
parameters {
  real drugV[D,K];
  vector[C] cellV[K];
  real<lower=0> noise; 
  matrix[D,D] meanSyn;
  real<lower=0> drugscale[K];
  real<lower=0> scalescale;
  real<lower=0> eta_sq;
  real<lower=0> inv_rho_sq;
  real<lower=0> sigma_sq;
}
transformed parameters {
      real<lower=0> rho_sq;
      rho_sq <- inv(inv_rho_sq);
}
model {
  matrix[C,C] Sigma;

  eta_sq ~ cauchy(0,5);
  inv_rho_sq ~ cauchy(0,5);
  sigma_sq ~ cauchy(0,5);


  for (i in 1:(C-1)) {
    for (j in (i+1):C) {
      Sigma[i,j] <- eta_sq * exp(-rho_sq * sqDist[i,j]);
      Sigma[j,i] <- Sigma[i,j];
    } 
  }
  for (k in 1:C)
    Sigma[k,k] <- eta_sq + sigma_sq;
  for (k in 1:K)
    cellV[k] ~ multi_normal(mu,Sigma);

  scalescale ~ cauchy(0,1);
  drugscale ~ cauchy(0,scalescale);
  for (d in 1:D)
    for (k in 1:K)
      drugV[d,k] ~ cauchy(0,drugscale[k]);
  
  noise ~ cauchy(0,5);
  // to_vector(meanSyn) ~ cauchy(0,5);
  for (n in 1:N) {
    real temp[K]; 
    for (k in 1:K)
      temp[k] <- drugV[drugA[n],k] * drugV[drugB[n],k] * cellV[k][cellLine[n] ]; 
    synergy[n] ~ normal( meanSyn[drugA[n],drugB[n]] + sum(temp), noise );
  }
}
generated quantities {
  real test_synergy[Ntest];
  for (n in 1:Ntest) {
    real temp[K]; 
    for (k in 1:K)
      temp[k] <- drugV[drugAtest[n],k] * drugV[drugBtest[n],k] * cellV[k][cellLinetest[n] ]; 
    test_synergy[n] <- meanSyn[drugA[n],drugB[n]] + sum(temp);
  }
}