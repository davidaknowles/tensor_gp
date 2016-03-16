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
  real y[N];
}
parameters {
  real drugV[D,K];
  matrix[C,K] cellV;
  real<lower=0> noise; 
  real meanSyn;
  real<lower=0> drugscale[K];
  real<lower=0> scalescale;
}
model {

  scalescale ~ cauchy(0,1);
  drugscale ~ cauchy(0,scalescale);
  for (d in 1:D)
    for (k in 1:K)
      drugV[d,k] ~ cauchy(0,drugscale[k]);
  
  noise ~ cauchy(0,5);

  to_vector( cellV ) ~ cauchy(0,1);

  for (n in 1:N) {
    real temp[K]; 
    for (k in 1:K)
      temp[k] <- drugV[drug[n],k] * cellV[cellLine[n],k ]; 
    y[n] ~ normal( meanSyn + sum(temp), noise );
  }
}
generated quantities {
  real test_synergy[Ntest];
  for (n in 1:Ntest) {
    real temp[K]; 
    for (k in 1:K)
      temp[k] <- drugV[drugtest[n],k] * cellV[cellLinetest[n],k]; 
    test_synergy[n] <- meanSyn + sum(temp);
  }
}