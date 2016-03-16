data {
  int<lower=0> K; 
  int<lower=0> K2; // hidden layers in NN
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
  vector[C] z[K];
  real<lower=0> noise; 
  real meanSyn;
  real<lower=0> drugscale[K];
  real<lower=0> scalescale;
  real<lower=0> eta_sq;
  real<lower=0> inv_rho_sq;
  real<lower=0> sigma_sq;
  matrix[K2,3*K] combBeta; // NN weights 
  vector[K2] combBeta2; 
  vector[K2] hiddenMu; 
}
transformed parameters {
  vector[C] cellV[K];
  matrix[C,C] Sigma;
  matrix[C,C] L;
  real<lower=0> rho_sq;
  rho_sq <- inv(inv_rho_sq);
  for (i in 1:(C-1)) {
    for (j in (i+1):C) {
      Sigma[i,j] <- eta_sq * exp(-rho_sq * sqDist[i,j]);
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
  eta_sq ~ cauchy(0,5);
  inv_rho_sq ~ cauchy(0,5);
  sigma_sq ~ cauchy(0,5);

  scalescale ~ cauchy(0,1);
  drugscale ~ cauchy(0,scalescale);
  for (d in 1:D)
    for (k in 1:K)
      drugV[d,k] ~ cauchy(0,drugscale[k]);
  
  noise ~ cauchy(0,5);
  
  to_vector( combBeta ) ~ cauchy(0,5);
  combBeta2 ~ cauchy(0,5);
  hiddenMu ~ cauchy(0,5);

  for (n in 1:N) {
    vector[3*K] temp; 
    vector[K2] hidden; 
    vector[K2] inputToTanh; 
    for (k in 1:K) {
      temp[k] <- drugV[drugA[n],k];
      temp[K+k] <- drugV[drugB[n],k];
      temp[2*K+k] <- cellV[k][cellLine[n] ];
    }
    inputToTanh <- combBeta * temp + hiddenMu; 
    for (k2 in 1:K2)
      hidden[k2] <- tanh( inputToTanh[k2] ); 
    synergy[n] ~ normal( meanSyn + dot_product(combBeta2 , hidden) , noise );
  }
}
generated quantities {
  real test_synergy[Ntest];
  for (n in 1:Ntest) {
    real temp[K]; 
    for (k in 1:K)
      temp[k] <- drugV[drugAtest[n],k] * drugV[drugBtest[n],k] * cellV[k][cellLinetest[n] ]; 
    test_synergy[n] <- meanSyn + sum(temp);
  }
}