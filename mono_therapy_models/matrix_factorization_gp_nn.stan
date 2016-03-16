data {
  int<lower=0> K; 
  int<lower=0> K2;
  int<lower=0> D; 
  int<lower=0> C;
  int<lower=0> N;
  int<lower=0,upper=C> drug[N]; 
  int<lower=0,upper=C> cellLine[N]; 
  int<lower=0> Ntest;
  int<lower=0,upper=C> drugtest[Ntest]; 
  int<lower=0,upper=C> cellLinetest[Ntest];
  matrix[C,C] sqDist; 
  vector[N] y;
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
  matrix[K2,2*K] combBeta; // NN weights 
  vector[K2] combBeta2; 
  vector[K2] hiddenMu; 
}
transformed parameters {
  vector[C] cellV[K];
  matrix[C,C] Sigma;
  matrix[C,C] L;
  real<lower=0> rho_sq;
  vector[2*K] x[N];
  vector[2*K] xtest[Ntest];
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
  for (n in 1:N) {
    for (k in 1:K) {
      x[n][k] <- drugV[drug[n],k];
      x[n][K+k] <- cellV[k][cellLine[n] ];
    }
  }
  for (n in 1:Ntest) {
    for (k in 1:K) {
      xtest[n][k] <- drugV[drugtest[n],k];
      xtest[n][K+k] <- cellV[k][cellLinetest[n] ];
    }
  }
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

  to_vector( combBeta ) ~ cauchy(0,1);
  combBeta2 ~ cauchy(0,1);
  hiddenMu ~ cauchy(0,1);

  for (n in 1:N) {
    vector[K2] hidden; 
    vector[K2] inputToTanh; 
    inputToTanh <- combBeta * x[n] + hiddenMu; 
    for (k2 in 1:K2)
      hidden[k2] <- tanh( inputToTanh[k2] ); 
    y[n] ~ normal( meanSyn + dot_product(combBeta2 , hidden) , noise );
  }
}
generated quantities {
  vector[Ntest] testMean;
  for (n in 1:Ntest) {
    vector[K2] hidden; 
    vector[K2] inputToTanh; 
    inputToTanh <- combBeta * xtest[n] + hiddenMu; 
    for (k2 in 1:K2)
      hidden[k2] <- tanh( inputToTanh[k2] ); 
    testMean[n] <- meanSyn + dot_product(combBeta2 , hidden) ;
  }
}